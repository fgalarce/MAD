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

#include<wave.hpp>

void Wave::initialize(Parameters parameters, const Geometry & geometry, Boundary & boundary){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "Wave: Initializing" << endl;
  par = parameters;
  geo = geometry;
  bd = boundary;
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

  mat(A, nbDofs, nbDofs); 
  mat(M, nbDofs, nbDofs); 

  fem.initialize(par, geo, bd, A);
  femMass.initialize(par, geo, bd, M);

}

void Wave::finalize(){
  MatDestroy(&M);
  VecDestroy(&u0);
  VecDestroy(&u00);
  VecDestroy(&dt_u);
  VecDestroy(&dt_u0);
  VecDestroy(&dtt_u);
  VecDestroy(&dtt_u0);
  KSPDestroy(&ksp);
}

Mat Wave::assembleLHS(){

  for (int partId = 0; partId < geo.elements().size(); partId++){
    for (vector<int> simplex : geo.elements()[partId]){

      fem.setSimplex(simplex);      
      femMass.copySimplex(fem);

      fem.grad_u_grad_v(par.youngModulus()[partId]/(2*(1.0+par.poissonRatio()[partId])));
      femMass.u_dot_v(par.density());

    }
  }

  code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERR(code);

  /* Rayleigh damping */
  double alpha = 0.47;
  double beta = 0.04;
  code = MatAXPY(A, beta, A, DIFFERENT_NONZERO_PATTERN); CHKERR(code);
  code = MatAXPY(A, alpha, M, DIFFERENT_NONZERO_PATTERN); CHKERR(code);

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
  if (m_world_rank == 0) cout << "Wave: Norm LHS = " << normLHS << endl;

  return K;
}

Vec Wave::assembleRHS(){
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
  if (m_world_rank == 0) cout << "Wave: Norm RHS = " << normRHS << endl;
  return b;
}

void Wave::update(Vec u, double time){
  if (par.timeIntegration() == "centeredFD" ){
    if (time > 0){
      code = VecZeroEntries(u00); CHKERR(code);
      code = VecAXPY(u00, 1.0, u0); CHKERR(code);
    }
    code = VecZeroEntries(u0); CHKERR(code);
    code = VecAXPY(u0, 1.0, u); CHKERR(code);
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

void Wave::setSolver(){
  if (m_world_rank == 0) cout << "KSP: Configuring" << endl;
  configureKSP(ksp, par.solver(), par.preconditioner(), par.monitorKSP(), par.use_solution_as_guess_KSP(), par.reuse_preconditioner());
}

void Wave::setLHS(Mat C){
  if (m_world_rank == 0) cout << "KSP: Seting operators" << endl;
  code = KSPSetOperators(ksp, C, C); CHKERR(code); 
  code = KSPSetUp(ksp); CHKERR(code);
}

void Wave::solve(Vec b, Vec u){
  if (m_world_rank == 0) cout << "Wave: Solving linear system" << endl;
  code = KSPSolve(ksp, b, u); CHKERR(code); 
  int its;
  double rnorm;
  code = KSPGetResidualNorm(ksp, &rnorm); CHKERR(code);
  code = KSPGetIterationNumber(ksp, &its); CHKERR(code);
  if (par.solver() == "gmres"){
    if (m_world_rank == 0) cout << "Wave: " << its << " iterations. Residual norm = " << rnorm << endl;
  }
}
