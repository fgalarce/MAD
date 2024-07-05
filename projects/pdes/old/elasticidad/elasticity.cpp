/*=============================================================================
  This file is part of the code MAD (or MAD)
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2021,
    
     Felipe Galarce at INRIA

  MDMA is free software; you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by the Free
  Software Foundation; either version 2.1 of the License, or (at your option)
  any later version.

  MDMA is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with MDMA. If not, see http://www.gnu.org/licenses/.
  =============================================================================*/

#include "assemble.hpp"

vector<double> inlet(vector<double> x, double t, Parameters par){
  vector<double> bc(2, 0.0);
  vector<double> center(2,0.0);
  double r = 4.0;
  bc[0] = - par.amplitude() * (dot(center - x, center - x)/(r*r) - 1.0);
//  bc[0] = par.amplitude();
  return bc;
}

vector<double> outlet(vector<double> x, double t, Parameters par){
  vector<double> bc(2, 0.0);
  return bc;
}

vector<double> wall(vector<double> x, double t, Parameters par){
  vector<double> bc(1, par.csf_pressure());
  return bc;
}

int main(int argc, char *argv[]){

  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  /* Parse data file */
  assert(argc == 2);
  string data_file = argv[1]; 
  Parameters par(data_file);
  par.print();

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);
    
  Geometry geo;
  geo.initialize(par, io);

  MasterElement fe; 
  fe.initialize(par, io.dimension());

  Boundary boundary;
  boundary.initialize(par, geo);

  int nbDofsPerNode = par.nbDofsPerNode()[0]; 
  int nbVertices = io.nbVertices();
  int nbDofs = nbVertices*nbDofsPerNode;

  /* local values */
  int nbNodesPerElement = geo.dimension()+1; 

  Mat A = mat(nbDofs, nbDofs); /* Static matrix */
//
//  assembleSystem(A, par, geo, fe);
//
//  /* LHS static matrix */
//  Mat C = mat(nbDofs, nbDofs);
//  code = MatDuplicate(A, MAT_COPY_VALUES, &C); CHKERRQ(code);
//
//  /* Boundary conditions */
////  boundary.Neumann(par.bcNeumann()[0], inlet, 0);
//  boundary.Dirichlet(par.inlet(), inlet, 0);
//  boundary.Dirichlet(par.outlet(), outlet, 0);
//
////  /* Uniform pressure at boundaries  */
////  for (int i : par.walls()){
////    boundary.Dirichlet(i, wall, 1);
////  }
////  boundary.Dirichlet(par.outlet(), wall, 1);
////  boundary.Dirichlet(par.inlet(), wall, 1);
//
//  /* Krylov solver */
//  KSP ksp;
//  configureKSP(ksp, par.solver(), par.preconditioner(), par.monitorKSP(), par.use_solution_as_guess_KSP(), par.reuse_preconditioner());
//  boundary.block(C);
//
//  if (world_rank == 0) cout << "Elasticity: Seting operators" << endl;
//  code = KSPSetOperators(ksp, C, C); CHKERR(code); 
//  code = KSPSetUp(ksp); CHKERRQ(code);
//
//  /* block and set up system of equations */ 
//  Vec u = zeros(nbDofs);
//  Vec b = zeros(nbDofs);
//  boundary.block(b);
//  if (world_rank == 0) cout << "Elasticity: Solving linear system" << endl;
//  code = KSPSolve(ksp, b, u); CHKERRQ(code); 
//  int its;
//  double rnorm;
//  code = KSPGetResidualNorm(ksp, &rnorm); CHKERRQ(code);
//  code = KSPGetIterationNumber(ksp, &its); CHKERRQ(code);
//
//  if (world_rank == 0) cout << "Elasticity: " << its << " iterations. Residual norm = " << rnorm << endl;
//  
//  io.writeState(u, 0.0);
//
//  double norm_sol = norm(u);
//  if (world_rank == 0) cout << "Elasticity: Norm solution = " << norm_sol << endl;
//
//  KSPDestroy(&ksp);
//  MatDestroy(&A);
//  VecDestroy(&u);
//  MatDestroy(&C);
  par.finalize();
  SlepcFinalize();
}
