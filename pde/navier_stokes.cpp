/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2023,
    
     Felipe Galarce at INRIA/WIAS/PUCV

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

#include<navier_stokes.hpp>

void NavierStokes::initialize(Parameters parameters, const Geometry & geometry, const Boundary & boundary){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "Navier Stokes: Initializing" << endl;
  par = parameters;
  geo = geometry;
  bd = boundary;
  calculus.initialize(par, geo, bd);
  fe.initialize(par, geo.dimension());
  feBD.initialize(par, geo.dimension());
  m_verbose = par.verbose();
  nbDofs = 0;
  assert(par.nbDofsPerNode().size() == 2);
  m_interpolator.initialize(par, geo);

  nbDofs = (par.nbDofsPerNode()[0] + par.nbDofsPerNode()[1])*geo.nbVertices;
  nbDofsPress = geo.nbVertices;
  nbDofsVel = par.nbDofsPerNode()[0]*geo.nbVertices;
  nbVertices = geo.nbVertices;
  nbDofsPerNode = par.nbDofsPerNode()[0]; 
  nbNodesPerElement = geo.dimension()+1; 

  u0 = zeros(nbDofs);
  VecAssemblyBegin(u0);
  VecAssemblyEnd(u0);

  mat(A, nbDofs, nbDofs); 
//  mat(B, nbDofs, nbDofs); 
  mat(M, nbDofs, nbDofs); 
//  mat(K, nbDofs, nbDofs); 
  mat(C, nbDofs, nbDofs); 
//  mat(C_rhs, nbDofs, nbDofs); 

  femStat.initialize(par, geo, bd, A);
//  femNonSym.initialize(par, geo, bd, B);
  femMass.initialize(par, geo, bd, M);
  femConv.initialize(par, geo, bd, C);
//  femRHS.initialize(par, geo, bd, C_rhs);
//  femStiff.initialize(par, geo, bd, K);
  code = MatGetOwnershipRange(A, &m, &n); CHKERR(code);

  if (!(par.power_law_n() == 1)){
    m_non_newtonian = true;
    if (m_world_rank == 0) cout << "Navier Stokes: The fluid is Non-Newtonian." << endl;
  }
  assembleInnerProduct();
}

void NavierStokes::finalize(){
  MatDestroy(&M);
  MatDestroy(&A);
  MatDestroy(&C);
  VecDestroy(&u0);
  KSPDestroy(&ksp);
}

Mat NavierStokes::assembleLHS_static(){

  for (int partId = 0; partId < geo.elements().size(); partId++){
    if (m_world_rank == 0) cout << "Navier Stokes: Assembling discretization matrix. Part " << partId << endl;
    for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on elements */

      vector<int> simplex = geo.elements()[partId][feId];

      femStat.setSimplex(simplex);
      femMass.copySimplex(femStat); 
//      femStiff.copySimplex(femStat);
//      femNonSym.copySimplex(femStat);

      if (feId % (geo.elements()[partId].size() / m_verbose) == 0){
        if (m_world_rank == 0) cout << "  Elementary matrix for element: " << feId << "/" << geo.elements()[partId].size() - 1 << endl;
      }

      /* u \cdot v */
      femMass.u_dot_v(par.density(), 0, 0);

      /* mu nabla u : nabla v */
      if (m_non_newtonian == false){ femStat.epsilon_u_epsilon_v(2.0*par.viscosity(), 0, 0); }
//      if (m_non_newtonian == false){ femStat.grad_u_grad_v(2.0*par.viscosity(), 0, 0); }

      /* -p div v */
      femStat.u_div_v(-1.0, 0, 1);

      /* div u q */
      femStat.div_u_v(1.0, 1, 0);

      /* SUPG  */
      double tau = femStat.feSize()*femStat.feSize();
      femStat.grad_u_grad_v(tau, 1, 1);

//      /* \nabla u \cdot \nabla v*/
//      femStiff.grad_u_grad_v(femStat.feSize()*femStat.feSize(), 0, 0);
//      femStat.grad_u_grad_v(tau*2*par.density()*par.viscosity()/par.timeStep(), 0, 0);
//      femNonSym.u__dot__grad_q(par.density()/par.timeStep()*tau, 1, 0);

    }
  }
  code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERR(code);
//  code = MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY); CHKERR(code); 
//  code = MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY); CHKERR(code);
//  code = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY); CHKERR(code); 
//  code = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY); CHKERR(code);

  /* LHS static matrix */
  Mat S = mat(nbDofs, nbDofs);
  code = MatDuplicate(M, MAT_COPY_VALUES, &S); CHKERR(code);
  code = MatScale(S, 1.0/par.timeStep());
  code = MatAXPY(S, 1.0, A, DIFFERENT_NONZERO_PATTERN); CHKERR(code);

  double normLHS = norm(S);
  if (m_world_rank == 0) cout << "Navier Stokes: Norm LHS = " << normLHS << endl;
  
  return S;
}

Vec NavierStokes::assembleRHS(Vec u0_custom){
  Vec b = zeros(nbDofs);
  Vec u0_seq;
  if (u0_custom == NULL){
    code = MatMult(M, u0, b); CHKERR(code);
    u0_seq = getSequential(u0);
  } else {
    code = MatMult(M, u0_custom, b); CHKERR(code);
    u0_seq = getSequential(u0_custom);
  }
  VecScale(b, 1.0/par.timeStep());

//  code = MatZeroEntries(C_rhs); CHKERR(code);
//
//  for (int partId = 0; partId < geo.elements().size(); partId++){
//    if (m_world_rank == 0) cout << "Navier Stokes: Assembling RHS. Part " << partId << endl;
//    for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on elements */
//
//      vector<int> simplex = geo.elements()[partId][feId];
//      vector<double> u_element = geo.getNodalValues(u0_seq, simplex, par.nbDofsPerNode()[0]);
//
//      femRHS.setSimplex(simplex);
//      double tau = femRHS.feSize()*femRHS.feSize();
//      femRHS.u__dot__a_grad_v(par.density()/par.timeStep()*tau, u_element, 0, 0);
//
//    }
//  }
//  MatAssemblyBegin(C_rhs, MAT_FINAL_ASSEMBLY);
//  MatAssemblyBegin(C_rhs, MAT_FINAL_ASSEMBLY);
//
//  /* SUPG stabilization RHS */
//  Vec b_supg = zeros(nbDofs);
//  if (u0_custom == NULL){
//    Vec b_supg1 = zeros(nbDofs);
//    Vec b_supg2 = zeros(nbDofs);
//    Vec b_supg3 = zeros(nbDofs);
//    code = MatMult(K, u0, b_supg1); CHKERR(code);
//    code = MatMult(B, u0, b_supg2); CHKERR(code);
//    code = MatMult(C_rhs, u0, b_supg3); CHKERR(code);
//    code = VecAXPY(b_supg, -2.0*par.viscosity()*par.density(), b_supg1); CHKERR(code);
//    code = VecAXPY(b_supg, -1.0, b_supg2); CHKERR(code);
//    code = VecAXPY(b_supg,  1.0, b_supg3); CHKERR(code);
//  } else {
//    Vec b_supg1 = zeros(nbDofs);
//    Vec b_supg2 = zeros(nbDofs);
//    Vec b_supg3 = zeros(nbDofs);
//    code = MatMult(K, u0_custom, b_supg1); CHKERR(code);
//    code = MatMult(B, u0_custom, b_supg2); CHKERR(code);
//    code = MatMult(C_rhs, u0_custom, b_supg3); CHKERR(code);
//    code = VecAXPY(b_supg, -2.0*par.viscosity()*par.density(), b_supg1); CHKERR(code);
//    code = VecAXPY(b_supg, -1.0, b_supg2); CHKERR(code);
//    code = VecAXPY(b_supg,  1.0, b_supg3); CHKERR(code);
//  }
//  if (!m_non_newtonian){
//    code = VecAXPY(b, -2.0*par.viscosity()*par.density(), b_supg); CHKERR(code);
//  } else {
//    exit(1);
//  }
//  VecDestroy(&b_supg);
  double normRHS = norm(b);
  if (m_world_rank == 0) cout << "Navier Stokes: Norm RHS = " << normRHS << endl;
  return b;
}

void NavierStokes::update(Vec u, double time){
  code = VecZeroEntries(u0); CHKERR(code);
  code = VecAXPY(u0, 1.0, u); CHKERR(code);
}

void NavierStokes::computeFlows(Vec u){
  double incompressibility_failure = 0.0;
  for (int i = 0; i < geo.bdLabels().size(); i++){
    double Q_bd = bd.flow(u, geo.bdLabels()[i]);
    if (m_world_rank == 0) cout << "Navier-Stokes: flow at boundary label " << geo.bdLabels()[i] << " : " << Q_bd << endl; 
    incompressibility_failure = incompressibility_failure + Q_bd;
  }
  if (m_world_rank == 0) cout << "Navier-Stokes: incompressibility failure = " << incompressibility_failure << endl; 
}

void NavierStokes::setSolver(){
  if (m_world_rank == 0) cout << "KSP: Configuring" << endl;
  configureKSP(ksp, par.solver(), par.preconditioner(), par.monitorKSP(), par.use_solution_as_guess_KSP(), par.reuse_preconditioner());
}

void NavierStokes::setLHS(Mat A, Mat P){
  if (m_world_rank == 0) cout << "KSP: Seting operators" << endl;
  if (P == NULL){
    code = KSPSetOperators(ksp, A, A); CHKERR(code); 
  } else {
    code = KSPSetOperators(ksp, A, P); CHKERR(code); 
  }
  code = KSPSetUp(ksp); CHKERR(code);
}

void NavierStokes::solve(Vec b, Vec u){
  if (m_world_rank == 0) cout << "NavierStokes: Solving linear system" << endl;
  code = KSPSolve(ksp, b, u); CHKERR(code); 
  
  int its;
  double rnorm;
  code = KSPGetResidualNorm(ksp, &rnorm); CHKERR(code);
  code = KSPGetIterationNumber(ksp, &its); CHKERR(code);
  if (par.solver() == "gmres"){
    if (m_world_rank == 0) cout << "NavierStokes: " << its << " iterations. Residual norm = " << rnorm << endl;
  }
}

Mat NavierStokes::assembleLHS(Mat LHS_static, Vec u0_custom){

  Vec u0_seq;
  double u0_norm;
  if (u0_custom == NULL){
    u0_seq = getSequential(u0);
    u0_norm = sqrt(innerProduct(calculus.split(u0)[0],calculus.split(u0)[0]));
  } else {
    u0_seq = getSequential(u0_custom);
    u0_norm = sqrt(innerProduct(calculus.split(u0)[0],calculus.split(u0)[0]));
  }
  code = MatZeroEntries(C); CHKERR(code);

  for (int partId = 0; partId < geo.elements().size(); partId++){
    if (m_world_rank == 0) cout << "Navier Stokes: Assembling convection matrix. Part " << partId << endl;
    for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on elements */

      vector<int> simplex = geo.elements()[partId][feId];
      femConv.setSimplex(simplex);
  
      /* Non Newtonian apparent viscosity */
      if (m_non_newtonian){ femConv.grad_u_grad_v(m_viscosity[feId], 0, 0); }

      /* Semi-implicit convection */
      vector<double> u_element = geo.getNodalValues(u0_seq, simplex, par.nbDofsPerNode()[0]);
      femConv.a_grad_u_dot_v(par.density(), u_element, 0, 0);

      /* SUPG stabilization */
//      double tau = femConv.feSize()*femConv.feSize();
//      femConv.a_grad_u__dot__a_grad_v(par.density()*par.density()*tau, u_element, 0, 0);
//      femConv.a_grad_v__dot__grad_p(par.density()*tau, u_element, 0, 1);
//      femConv.a_grad_u__dot__grad_q(par.density()*tau, u_element, 1, 0);
//      femConv.u__dot__a_grad_v(par.density()/par.timeStep()*tau, u_element, 0, 0);
//      femConv.u__dot__grad_q(par.density()/par.timeStep()*tau, 0, 0);
    }
  }
  VecDestroy(&u0_seq);
  code = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY); CHKERR(code);

  /* assemble full LHS matrix */ 
  Mat S = mat(nbDofs, nbDofs);
  code = MatDuplicate(LHS_static, MAT_COPY_VALUES, &S); CHKERR(code);
  code = MatAXPY(S, 1.0, C, DIFFERENT_NONZERO_PATTERN); CHKERR(code);
  return S;
}

/* compute \grad u \cdot n */
Vec NavierStokes::shearStress(Vec u, int bdLabel){
  if (m_world_rank == 0) cout << "Navier Stokes: computing shear stress at boundary label " << bdLabel << "." << endl;
  return calculus.gradientOnBoundary(u, bdLabel);
}

vector<double> NavierStokes::getViscosity(Vec u0_custom){

  Vec u0_seq;
  if (u0_custom == NULL){
    u0_seq = getSequential(u0);
  } else {
    u0_seq = getSequential(u0_custom);
  }

  m_viscosity.resize(geo.elements()[0].size());
  vector<Vec> u0_seq_percomponent(par.nbDofsPerNode()[0]);

  for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
    u0_seq_percomponent[comp] = getSequential(calculus.decompose_vector_field(u0_seq)[comp]);
  }

  for (int partId = 0; partId < geo.elements().size(); partId++){
    if (m_world_rank == 0) cout << "Navier Stokes: Computing apparent viscosity. Part " << partId << endl;
    for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on elements */

      vector<int> simplex = geo.elements()[partId][feId];
      /* get finite element coordinates */
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int i = 0; i < nbNodesPerElement; i++){
        coordinates[i] = geo.coordinates()[simplex[i]];
      }
      fe.setCoordinates(coordinates);
      fe.computeSize();

      if (feId % (geo.elements()[partId].size() / par.verbose()) == 0){
        if (m_world_rank == 0) cout << "  Elementary vector for simplex: " << feId << "/" << geo.elements()[partId].size() - 1 << endl;
      }
   
      vector<vector<double>> u_element(par.nbDofsPerNode()[0]);
      vector<vector<double>> grad_u(par.nbDofsPerNode()[0]);
      for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
        u_element[comp] = geo.getNodalValues(u0_seq_percomponent[comp], simplex, 1);;
        grad_u[comp].resize(par.nbDofsPerNode()[0]);
      }
      for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
        for (int nodeId = 0; nodeId < nbNodesPerElement; nodeId++){
          grad_u[comp] = grad_u[comp] + u_element[comp][nodeId] * transpose(fe.jacobianInverse()) * fe.grad_phi(nodeId);
        }
      }
      grad_u = transpose(grad_u);

      double eta;
      double gamma_dot = sqrt(trace(grad_u*transpose(grad_u)));
      if (gamma_dot != 0){
        eta = par.power_law_m() * pow(gamma_dot, par.power_law_n() - 1.0);
      } 
      if (eta < par.power_law_m()/10000){
        eta = par.power_law_m()/10000;
      } else if (eta > 10000*par.power_law_m()){
        eta = par.power_law_m()*10000;
      }
      m_viscosity[feId] = eta;
    }
  }
  m_p1Viscosity = m_interpolator.interpolateP0_P1(m_viscosity);

  return m_viscosity;
}

double NavierStokes::innerProduct(Vec u, Vec v){
  int u_size; double ip;
  code = VecGetSize(u, &u_size); CHKERR(code);
  Vec Mv = vec(u_size);
  if (u_size == par.nbDofsPerNode()[0]*geo.nbVertices){
    code = MatMult(m_ip[0], v, Mv); CHKERR(code);
  } else if (u_size == par.nbDofsPerNode()[1]*geo.nbVertices){
    code = MatMult(m_ip[1], v, Mv); CHKERR(code);
//  } else if (u_size == (par.nbDofsPerNode()[0] + par.nbDofsPerNode()[1])*geo.nbVertices){
//    code = MatMult(m_ip[0], v, Mv); CHKERR(code);
//    Vec Mv2 = vec(u_size);
//    code = MatMult(m_ip[1], v, Mv2); CHKERR(code);
//    code = VecAXPY(Mv, 1.0, Mv2); CHKERR(code);
//    code = VecDestroy(&Mv2); CHKERR(code);
  } else {
    errorMessage("innerProduct", "Wrong vector size");
  }
  code = VecDot(Mv, u, &ip); CHKERR(code);
  code = VecDestroy(&Mv); CHKERR(code);

  return ip;
}

void NavierStokes::assembleInnerProduct(){

  m_ip.resize(2);
  mat(m_ip[0], nbDofsVel, nbDofsVel);
  mat(m_ip[1], nbDofsPress, nbDofsPress);

  for (int partId = 0; partId < geo.elements().size(); partId++){
    if (m_world_rank == 0) cout << "Navier Stokes: Assembling inner product matrix. Part " << partId << endl;
    for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on elements */

      vector<int> simplex = geo.elements()[partId][feId];
      /* get finite element coordinates */
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int i = 0; i < nbNodesPerElement; i++){
        coordinates[i] = geo.coordinates()[simplex[i]];
      }

      fe.setCoordinates(coordinates);
      fe.computeSize();

      if (feId % (geo.elements()[partId].size() / m_verbose) == 0){
        if (m_world_rank == 0) cout << "  Elementary matrix for element: " << feId << "/" << geo.elements()[partId].size() - 1 << endl;
      }
      /* Assemble elementary matrices */
      for (int j = 0; j < nbNodesPerElement; j++){
        for (int i = 0; i < nbNodesPerElement; i++){
          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
            if (par.nbDofsPerNode()[0]*simplex[i]+comp >= m && par.nbDofsPerNode()[0]*simplex[i]+comp < n){
              double Kij = fe.stiffness(i,j);
              double Mij = fe.mass(i,j);
              code = MatSetValue(m_ip[0], par.nbDofsPerNode()[0]*simplex[i]+comp, par.nbDofsPerNode()[0]*simplex[j]+comp, Kij+Mij, ADD_VALUES); CHKERR(code);
            }
          }
          if (nbDofsVel + simplex[i] >= m && nbDofsVel + simplex[i] < n){
            double Mij = fe.mass(i,j);
            code = MatSetValue(m_ip[1], simplex[i], simplex[j], Mij, ADD_VALUES); CHKERR(code); /* B^T */ 
          }
        }
      }
    }
  }
  code = MatAssemblyBegin(m_ip[0], MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(m_ip[0], MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyBegin(m_ip[1], MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(m_ip[1], MAT_FINAL_ASSEMBLY); CHKERR(code);
}
