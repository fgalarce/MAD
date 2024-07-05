/*=============================================================================
  This file is part of the code MAD
  Copyright (C) 2017-2024,
    
     Felipe Galarce

  MAD is free software; you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by the Free
  Software Foundation; either version 2.1 of the License, or (at your option)
  any later version.

  MAD is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with MAD. If not, see http://www.gnu.org/licenses/.
  =============================================================================*/

#include<hammer.hpp>

void Hammer::initialize(Parameters parameters, const Geometry & geometry){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "Hammer: Initializing" << endl;
  par = parameters;
  geo = geometry;
  fe.initialize(par, geo.dimension());
  m_verbose = par.verbose();
  nbDofs = 0;
  assert(par.nbDofsPerNode().size() == 2);

  nbVertices = geo.nbVertices;
  nbDofs = (par.nbDofsPerNode()[0] + par.nbDofsPerNode()[1])*nbVertices;
  nbDofsPress = geo.nbVertices;
  nbDofsVel = par.nbDofsPerNode()[0]*geo.nbVertices;
  nbDofsPerNode = par.nbDofsPerNode()[0]; 

  /* Adimensionalisation */
  m_Pr = par.density()*9.81*par.heigth();
  m_tau0 = par.bingham() * m_Pr * par.diameter() / 4.0 / par.length();
  m_Ur = m_Pr*par.diameter()*par.diameter()/32.0/par.length()/par.viscosity();
  double ReRef = par.density() * par.diameter() * m_Ur / par.viscosity();
  double MaRef = m_Ur / par.wave_velocity()[0];
  double aspect_ratio = par.diameter()/par.length();
  m_lambda_1 = 1.0/32*ReRef*aspect_ratio;
  m_lambda_2 = 1.0/32*ReRef/MaRef*aspect_ratio;

  double He = par.diameter() * par.diameter() * m_tau0 * par.density() / par.viscosity() / par.viscosity();
  m_fReRef = 16 + (10.67 + 0.1414*pow(He/ReRef, 1.143)) / (1.0 + 0.0149*pow(He/ReRef, 1.16) ) * (He/ReRef);

  if (m_world_rank == 0) cout << "u*(t=0.0) = 16 / f Re_ref = " <<16 /  m_fReRef << endl;
  if (m_world_rank == 0) cout << "Adimensionalisation: Ur = " << m_Ur << " m/s " << ", Pr = " << m_Pr << " Pa, tau_y = " << m_tau0 << " Pa, He = " << He << endl; 

  /* Set initial conditions */
  u0 = zeros(nbDofs);
  VecSet(u0, 16.0/m_fReRef);
  double feSize = geo.coordinates()[1][0] - geo.coordinates()[0][0];
  for (int i = nbVertices; i < 2*nbVertices; i++){
    code = VecSetValue(u0, i, (nbVertices-i)*feSize, INSERT_VALUES); CHKERR(code);
  }
  VecAssemblyBegin(u0);
  VecAssemblyEnd(u0);

  mat(C, nbDofs, nbDofs); 
}

void Hammer::compute_fRe(double velocity, int feId){
  double Re = par.density() * fabs(velocity) * m_Ur * par.diameter() / par.viscosity();
  double He = par.diameter() * par.diameter() * m_tau0 * par.density() / par.viscosity() / par.viscosity();
  double frictionRe = 16.0;
  if (Re < 4000 && Re > 1e-1){
    /* Buckingham–Reiner */
    frictionRe = 16 + (10.67 + 0.1414*pow(He/Re, 1.143)) / (1.0 + 0.0149*pow(He/Re, 1.16) ) * (He/Re);
  } else if (Re > 4000) {
    /* Darby and Melson */
    double a = -1.47 * (1.0 + 0.146 * exp(-2.9*1e-5*He));
    frictionRe = Re * pow(10, a) / pow(Re, 0.193);
  }
  fRe[feId] = frictionRe;
  Reynolds[feId] = Re;
}

void Hammer::finalize(){
  MatDestroy(&M);
  MatDestroy(&Mp);
  MatDestroy(&A);
  VecDestroy(&u0);
  KSPDestroy(&ksp);
}

Mat Hammer::assembleLHS(Mat LHS_static, Vec u0_custom){

  Vec u0_seq;
  if (u0_custom == NULL){
    u0_seq = getSequential(u0);
  } else {
    u0_seq = getSequential(u0_custom);
  }
  code = MatZeroEntries(C); CHKERR(code);
  
  /* local values */
  int nbNodesPerElement = geo.dimension()+1; 

  int m,n;
  code = MatGetOwnershipRange(C, &m, &n); CHKERR(code);
  code = MatZeroEntries(C); CHKERR(code);

  double ReRef = par.density() * par.diameter() * m_Ur / par.viscosity();
  for (int partId = 0; partId < geo.elements().size(); partId++){

    if (m_world_rank == 0) cout << "Hammer: Assembling discretization matrix. Part " << partId << endl;
    if (m_world_rank == 0) cout << "Adimensionalisation: Ur = " << m_Ur << " m/s " << ", Pr = " << m_Pr << " Pa, tau_y = " << m_tau0 << ", fReRef = " << m_fReRef << endl; 
    if (m_world_rank == 0) cout << "lambda_1 = " << m_lambda_1 << endl; 
    if (m_world_rank == 0) cout << "lambda_2 = " << m_lambda_2 << endl; 
    if (m_world_rank == 0) cout << "Average Reynolds in space = " << average(Reynolds) << endl; 
    fRe.resize(geo.elements()[partId].size(), 0.0);
    Reynolds.resize(geo.elements()[partId].size(), 0.0);

    for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on finiteElement */
      
      /* get finite element coordinates */
      vector<int> simplex = geo.elements()[partId][feId];
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int nodeId = 0; nodeId < nbNodesPerElement; nodeId++){
        coordinates[nodeId] = geo.coordinates()[simplex[nodeId]];
      }
      fe.setCoordinates(coordinates);
      fe.computeSize();
      /* Assemble elementary matrices */
      vector<double>  u_element = geo.getNodalValues(u0_seq, simplex, 1);
      double mean_velocity = (u_element[0] + u_element[1])/2.0;
      compute_fRe(mean_velocity, feId);
      for (int i = 0; i < nbNodesPerElement; i++){
        for (int j = 0; j < nbNodesPerElement; j++){
          if (simplex[i] >= m && simplex[i] < n){
            double Mij = fe.mass(i,j);
            code = MatSetValue(C, simplex[i], simplex[j], fRe[feId] / 16.0 * Mij, ADD_VALUES); CHKERR(code);

            /* Convection */
            if (par.convective_hammer()){
              double Cij = fe.phi_dphi_i_COMP_phi_j(j,i,0,u_element); 
              code = MatSetValue(C, simplex[i], simplex[j], m_lambda_1*Cij, ADD_VALUES); CHKERR(code);
            }
          }
        }
      }
    }
  }
  code = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY); CHKERR(code);

  double normLHS = norm(C);
  if (m_world_rank == 0) cout << "Hammer: Norm LHS time dependent = " << normLHS << endl;

  code = MatAXPY(C, 1.0, LHS_static, DIFFERENT_NONZERO_PATTERN); CHKERR(code);

  normLHS = norm(C);
  if (m_world_rank == 0) cout << "Hammer: Norm LHS = " << normLHS << endl;

  return C;
}

Mat Hammer::assembleLHS_static(){
  mat(A, nbDofs, nbDofs); 
  mat(M, nbDofs, nbDofs); 

  /* local values */
  int nbNodesPerElement = geo.dimension()+1; 
  int m,n;
  code = MatGetOwnershipRange(A, &m, &n); CHKERR(code);
  for (int partId = 0; partId < geo.elements().size(); partId++){
    if (m_world_rank == 0) cout << "Hammer: Assembling discretization matrix. Part " << partId << endl;
    for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on finiteElement */
      
      /* get finite element coordinates */
      vector<int> simplex = geo.elements()[partId][feId];
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int nodeId = 0; nodeId < nbNodesPerElement; nodeId++){
        coordinates[nodeId] = geo.coordinates()[simplex[nodeId]];
      }
      fe.setCoordinates(coordinates);
      fe.computeSize();

      if (feId % (geo.elements()[partId].size() / par.verbose()) == 0){
        if (m_world_rank == 0) cout << "    Elementary matrix for element: " << feId << "/" << geo.elements()[partId].size() - 1 << endl;
      }
      /* Assemble elementary matrices */
      for (int i = 0; i < nbNodesPerElement; i++){
        for (int j = 0; j < nbNodesPerElement; j++){
          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
            /* Momentum conservation */
            if (nbDofsPerNode*simplex[i]+comp >= m && nbDofsPerNode*simplex[i]+comp < n){
              double Mij = fe.mass(i,j);
              /* u \cdot v */
              code = MatSetValue(M, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp, (m_lambda_2/par.timeStep()) * Mij, ADD_VALUES); CHKERR(code);
              /* dp/dx v */
              double Bij = fe.mixed(j,i,comp);
              code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp, nbDofsVel + simplex[j], 1.0*Bij, ADD_VALUES); CHKERR(code); /* -B*/ 
            }
            /* Mass conservation */
            if (nbDofsVel + simplex[i] >= m && nbDofsVel + simplex[i] < n){
              double Bij = fe.mixed(j,i,comp);
              code = MatSetValue(A, nbDofsVel + simplex[i], nbDofsPerNode*simplex[j]+comp, m_lambda_2*Bij, ADD_VALUES); CHKERR(code); /* B^T */ 
            }
          }
          if (nbDofsVel + simplex[i] >= m && nbDofsVel + simplex[i] < n){
            double Mij = fe.mass(i,j);
            code = MatSetValue(M, nbDofsVel + simplex[i], nbDofsVel + simplex[j], 1.0/par.timeStep() * Mij, ADD_VALUES); CHKERR(code);
          }
        }
      }
    }
  }

  code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERR(code);

  /* LHS static matrix */
  code = MatAXPY(A, 1.0, M, DIFFERENT_NONZERO_PATTERN); CHKERR(code);

  double normLHS = norm(A);
  if (m_world_rank == 0) cout << "Hammer: Norm LHS = " << normLHS << endl;

  return A;
}

Vec Hammer::assembleRHS(Vec u0_custom){
  Vec b = zeros(nbDofs);
  if (u0_custom == NULL){
    code = MatMult(M, u0, b); CHKERR(code);
  } else {
    code = MatMult(M, u0_custom, b); CHKERR(code);
  }
  double normRHS  = norm(b);
  if (m_world_rank == 0) cout << "Hammer: Norm RHS = " << normRHS << endl;
  return b;
}

void Hammer::update(Vec u, double time){
  code = VecZeroEntries(u0); CHKERR(code);
  code = VecAXPY(u0, 1.0, u); CHKERR(code);
}

void Hammer::setSolver(){
  if (m_world_rank == 0) cout << "KSP: Configuring" << endl;
  configureKSP(ksp, par);
}

void Hammer::setLHS(Mat C){
  if (m_world_rank == 0) cout << "KSP: Seting operators" << endl;
  code = KSPSetOperators(ksp, C, C); CHKERR(code); 
  code = KSPSetUp(ksp); CHKERR(code);
}

void Hammer::solve(Vec b, Vec u){
  if (m_world_rank == 0) cout << "Hammer: Solving linear system" << endl;
  code = KSPSolve(ksp, b, u); CHKERR(code); 
  int its;
  double rnorm;
  code = KSPGetResidualNorm(ksp, &rnorm); CHKERR(code);
  code = KSPGetIterationNumber(ksp, &its); CHKERR(code);
  if (par.solver() == "gmres"){
    if (m_world_rank == 0) cout << "Hammer: " << its << " iterations. Residual norm = " << rnorm << endl;
  }
}
