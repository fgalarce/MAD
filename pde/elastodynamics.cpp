/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2021,
    
     Felipe Galarce at INRIA

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

#include<elastodynamics.hpp>

void Elastodynamics::initialize(Parameters parameters, const Geometry & geometry){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "Elastodynamics: Initializing" << endl;
  par = parameters;
  geo = geometry;
  fe.initialize(par, geo.dimension());
  m_verbose = par.verbose();
  nbDofs = 0;
  assert(par.nbDofsPerNode().size() == 1);

  nbDofsPerNode = par.nbDofsPerNode()[0]; 
  nbVertices = geo.nbVertices;
  nbDofs = nbVertices * nbDofsPerNode;

  assert(par.timeIntegration() == "newmark" || par.timeIntegration() == "centeredFD");

  u0 = zeros(nbDofs);
  if (par.timeIntegration() == "centeredFD"){
    u00 = zeros(nbDofs);
  } else if (par.timeIntegration() == "newmark"){
    delta = 1.0/2.0;
    alpha = 1.0/4.0;
    dt_u = zeros(nbDofs);
    dt_u0 = zeros(nbDofs);
    dtt_u = zeros(nbDofs);
    dtt_u0 = zeros(nbDofs);
    new_a0 = 1.0 / alpha / par.timeStep() / par.timeStep();
    new_a2 = 1.0 / alpha / par.timeStep();
    new_a3 = 1.0 / 2.0 / alpha - 1.0;
    new_a6 = par.timeStep() * (1.0 - delta);
    new_a7 = delta * par.timeStep();
  }
}

void Elastodynamics::finalize(){
  MatDestroy(&M);
  VecDestroy(&u0);
  VecDestroy(&u00);
  VecDestroy(&dt_u);
  VecDestroy(&dt_u0);
  VecDestroy(&dtt_u);
  VecDestroy(&dtt_u0);
  KSPDestroy(&ksp);
}

Mat Elastodynamics::assembleLHS(){

  mat(A, nbDofs, nbDofs); 
  mat(M, nbDofs, nbDofs); 

  /* local values */
  int nbNodesPerElement = geo.dimension()+1; 

  int m,n;
  code = MatGetOwnershipRange(A, &m, &n); CHKERR(code);

  for (int partId = 0; partId < geo.elements().size(); partId++){

    if (m_world_rank == 0) cout << "Elastodynamics: Assembling discretization matrix. Part " << partId << endl;
      /* Compute Lame coefficients */
    double lambda = par.youngModulus()[partId]/(2.0*(1.0+par.poissonRatio()[partId]));
    double mu = par.youngModulus()[partId] * par.poissonRatio()[partId] / ((1.0 + par.poissonRatio()[partId])*(1.0 - 2.0*par.poissonRatio()[partId]));
    if (m_world_rank == 0) cout << "Elastodynamics: Lame parameters. mu = " << mu << ". lambda = " << lambda << ". " << endl;

    MasterElement fe; 
    fe.initialize(par, geo.dimension());

    for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on finiteElement */
      
      /* get finite element coordinates */
      vector<int> simplex = geo.elements()[partId][feId];
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int nodeId = 0; nodeId < nbNodesPerElement; nodeId++){
        coordinates[nodeId] = geo.coordinates()[simplex[nodeId]];
      }
      fe.setCoordinates(coordinates);

      if (feId % (geo.elements()[partId].size() / par.verbose()) == 0){
        if (m_world_rank == 0) cout << "    Elementary matrix for element: " << feId << "/" << geo.elements()[partId].size() - 1 << endl;
        if (m_world_rank == 0) cout << "    detJac: " << fe.detJacobian() << endl;
      }
      /* Assemble elementary matrices */
      for (int i = 0; i < nbNodesPerElement; i++){
        for (int j = 0; j < nbNodesPerElement; j++){

          /* u \cdot v */
          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
            if (nbDofsPerNode*simplex[i]+comp >= m && nbDofsPerNode*simplex[i]+comp < n){
              double Mij = par.density() * fe.mass(i,j);
              code = MatSetValue(M, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp, Mij, ADD_VALUES); CHKERR(code);
            }
          }

          /* \epsilon(u) : \epsilon(v)  */
          for (int comp = 0; comp < nbDofsPerNode; comp++){
            if (nbDofsPerNode*simplex[i]+comp >= m && nbDofsPerNode*simplex[i]+comp < n){
              double Kij = 1.0*fe.stiffness_symmetric(i,j); /* diagonal components */
              code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp, mu*Kij, ADD_VALUES); CHKERR(code);
            }
          }
          for (int comp_v = 0; comp_v < nbDofsPerNode; comp_v++){
            if (nbDofsPerNode*simplex[i]+comp_v>= m && nbDofsPerNode*simplex[i]+comp_v < n){
              for (int comp_u = 0; comp_u < nbDofsPerNode; comp_u++){
                double Kdev = 1.0/2.0 * fe.dphi_dx(i, comp_u) * fe.dphi_dx(j, comp_v) * fe.detJacobian()*fe.weightsSum();
                code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp_v, nbDofsPerNode*simplex[j]+comp_u, mu*Kdev, ADD_VALUES); CHKERR(code);
              }
            }
          }

          /* div(u) div(v) */
          for (int comp_v = 0; comp_v < nbDofsPerNode; comp_v++){
            /* deviatoric components */
            if (nbDofsPerNode*simplex[i]+comp_v>= m && nbDofsPerNode*simplex[i]+comp_v < n){
              for (int comp_u = 0; comp_u < nbDofsPerNode; comp_u++){
                double div_div_ij = fe.dphi_dx(i, comp_v) * fe.dphi_dx(j, comp_u) * fe.detJacobian()*fe.weightsSum();
                code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp_v, nbDofsPerNode*simplex[j]+comp_u,  lambda*div_div_ij, ADD_VALUES); CHKERR(code);
              }
            }
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
  Mat K = mat(nbDofs, nbDofs);
  code = MatDuplicate(A, MAT_COPY_VALUES, &K); CHKERR(code);
  if (par.timeIntegration() == "centeredFD"){
    code = MatScale(M, 1.0/(par.timeStep()*par.timeStep())); CHKERR(code);
    code = MatAXPY(K, 1.0, M, DIFFERENT_NONZERO_PATTERN); CHKERR(code);
  } else if (par.timeIntegration() == "newmark") {
    code = MatAXPY(K, new_a0, M, DIFFERENT_NONZERO_PATTERN); CHKERR(code);
  }

  double normLHS = norm(K);
  if (m_world_rank == 0) cout << "Elastodynamics: Norm LHS = " << normLHS << endl;

  return K;
}

Vec Elastodynamics::assembleRHS(){
  Vec b = zeros(nbDofs);
  if (par.timeIntegration() == "centeredFD"){
    code = VecScale(u0, 2.0); CHKERR(code);
    code = VecAXPY(u0, -1.0, u00); CHKERR(code);
    code = MatMult(M, u0, b); CHKERR(code);
  } else if (par.timeIntegration() == "newmark") {
    Vec b0 = zeros(nbDofs);
    code = VecAXPY(b0, new_a0, u0); CHKERR(code);
    code = VecAXPY(b0, new_a2, dt_u); CHKERR(code);
    code = VecAXPY(b0, new_a3, dtt_u); CHKERR(code);
    code = MatMult(M, b0, b); CHKERR(code);
  }
  double normRHS  = norm(b);
  if (m_world_rank == 0) cout << "Elastodynamics: Norm RHS = " << normRHS << endl;
  return b;
}

void Elastodynamics::update(Vec u, double time){
  if (par.timeIntegration() == "centeredFD" ){
    code = VecZeroEntries(u0); CHKERR(code);
    code = VecAXPY(u0, 1.0, u); CHKERR(code);
    if (time > 0){
      code = VecZeroEntries(u00); CHKERR(code);
      code = VecAXPY(u00, 1.0, u0); CHKERR(code);
    }
  } else if (par.timeIntegration() == "newmark") {
    code = VecZeroEntries(dtt_u); CHKERR(code);
    /* re compute velocity */
    code = VecAXPY(dtt_u,  new_a0, u); CHKERR(code);
    code = VecAXPY(dtt_u, -new_a0, u0); CHKERR(code);
    code = VecAXPY(dtt_u, -new_a2, dt_u0); CHKERR(code);
    code = VecAXPY(dtt_u, -new_a3, dtt_u0); CHKERR(code);
    /* re compute acceleration */
    code = VecAXPY(dt_u, 1.0, dt_u0); CHKERR(code);
    code = VecAXPY(dt_u, new_a6, dtt_u0); CHKERR(code);
    code = VecAXPY(dt_u, new_a7, dtt_u); CHKERR(code);
    /* update states */
    code = VecZeroEntries(u0); CHKERR(code);
    code = VecAXPY(u0, 1.0, u); CHKERR(code);
    code = VecZeroEntries(dtt_u0); CHKERR(code);
    code = VecAXPY(dtt_u0, 1.0, dtt_u); CHKERR(code);
    code = VecZeroEntries(dtt_u); CHKERR(code);
    code = VecZeroEntries(dt_u0); CHKERR(code);
    code = VecAXPY(dt_u0, 1.0, dt_u); CHKERR(code);
    code = VecZeroEntries(dt_u); CHKERR(code);
  }
}

void Elastodynamics::setSolver(){
  if (m_world_rank == 0) cout << "KSP: Configuring" << endl;
  configureKSP(ksp, par);
}

void Elastodynamics::setLHS(Mat C){
  if (m_world_rank == 0) cout << "KSP: Seting operators" << endl;
  code = KSPSetOperators(ksp, C, C); CHKERR(code); 
  code = KSPSetUp(ksp); CHKERR(code);
}

void Elastodynamics::solve(Vec b, Vec u){
  if (m_world_rank == 0) cout << "Elastodynamics: Solving linear system" << endl;
  code = KSPSolve(ksp, b, u); CHKERR(code); 
  int its;
  double rnorm;
  code = KSPGetResidualNorm(ksp, &rnorm); CHKERR(code);
  code = KSPGetIterationNumber(ksp, &its); CHKERR(code);
  if (par.solver() == "gmres"){
    if (m_world_rank == 0) cout << "Elastodynamics: " << its << " iterations. Residual norm = " << rnorm << endl;
  }
}
