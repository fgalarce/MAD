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
//  bc[0] = - par.amplitude() * (dot(center - x, center - x)/(r*r) - 1.0);
  bc[0] = par.amplitude();
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
  int nbDofsSol = nbVertices * nbDofsPerNode;
  int nbDofsPress = nbVertices;
  int nbDofs = nbDofsSol + nbDofsPress;

  /* local values */
  int nbNodesPerElement = geo.dimension()+1; 

  /* initial condition */
  Vec u00 = zeros(nbDofs);
  Vec u0 = zeros(nbDofs);
  Vec u = zeros(nbDofs);

  //io.writeState(u0, 0.0);

  Mat M_sol_2 = mat(nbDofs, nbDofs); /* Mass matrix div */
  Mat A = mat(nbDofs, nbDofs); /* Static matrix */

  assembleSystem(A, M_sol_2, par, geo, fe);

  /* LHS static matrix */
  Mat C = mat(nbDofs, nbDofs);
  code = MatDuplicate(A, MAT_COPY_VALUES, &C); CHKERRQ(code);
  code = MatAXPY(C, 1.0, M_sol_2, DIFFERENT_NONZERO_PATTERN); CHKERRQ(code);

//  saveMat(C, "C.txt");

  /* Boundary conditions */
  boundary.Dirichlet(par.inlet(), inlet, 0);
  boundary.Dirichlet(par.outlet(), outlet, 0);

  /* Uniform pressure at boundaries  */
  for (int i : par.walls()){
    boundary.Dirichlet(i, wall, 1);
  }
  boundary.Dirichlet(par.outlet(), wall, 1);
  boundary.Dirichlet(par.inlet(), wall, 1);

  /* Krylov solver */
  KSP ksp;
  configureKSP(ksp, par.solver(), par.preconditioner(), par.monitorKSP(), par.use_solution_as_guess_KSP(), par.reuse_preconditioner());
  boundary.block(C);

  if (world_rank == 0) cout << "Poro: Seting operators" << endl;
  code = KSPSetOperators(ksp, C, C); CHKERR(code); 
  code = KSPSetUp(ksp); CHKERRQ(code);

  ofstream file_norm_solution(par.dirResults() + "/norm_solution.txt");
  
  for (double t : range(0, par.nbIterations()*par.timeStep(), par.timeStep())){

    if (world_rank == 0) cout << "\n- - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
    if (world_rank == 0) cout << "Poro: Solving poro-elasticity equation: it " <<  t/par.timeStep() << " t=" << t << endl;

    /* RHS Darcy */
    boundary.time(t); 
    Vec b_sol2 = zeros(nbDofs);
    code = MatMult(M_sol_2, u0, b_sol2); CHKERRQ(code);
    double Mu_sol2_norm = norm(b_sol2);
    if (world_rank == 0) cout << "Poro: Norm b_solid_2 = " << Mu_sol2_norm << endl;

    Vec b = zeros(nbDofs);
    code = VecAXPY(b, 1.0, b_sol2); CHKERRQ(code);

    /* Boundary conditions */
    boundary.Dirichlet(par.inlet(), inlet, 0);
    boundary.Dirichlet(par.outlet(), outlet, 0);

    for (int i : par.walls()){
      boundary.Dirichlet(i, wall, 1); /* confined , only axial movement allowed */
    }
    boundary.Dirichlet(par.outlet(), wall, 1);
    boundary.Dirichlet(par.inlet(), wall, 1);
   
    /* block and set up system of equations */ 
    boundary.block(b);
    if (world_rank == 0) cout << "Poro: Solving linear system" << endl;
    code = KSPSolve(ksp, b, u); CHKERRQ(code); 
    int its;
    double rnorm;
    code = KSPGetResidualNorm(ksp, &rnorm); CHKERRQ(code);
    code = KSPGetIterationNumber(ksp, &its); CHKERRQ(code);

    if (world_rank == 0) cout << "Poro: " << its << " iterations. Residual norm = " << rnorm << endl;
    
//    if ( t/par.timeStep() >= par.nbIterations() - 1){
      io.writeState(u, t + par.timeStep());
//    }
       
    /* Update solutions for centered finite differences */ 
    if (t > 0){
      VecZeroEntries(u00); 
      VecAXPY(u00, 1.0, u0);
    }
    code = VecZeroEntries(u0); CHKERRQ(code);
    code = VecAXPY(u0, 1.0, u); CHKERRQ(code);

    double norm_sol = norm(u);
    if (world_rank == 0) cout << "Poro: Norm solution = " << norm_sol << endl;
    file_norm_solution << norm_sol << endl; 

  }
  file_norm_solution.close();

  KSPDestroy(&ksp);
  MatDestroy(&A);
  VecDestroy(&u);
  VecDestroy(&u0);
  VecDestroy(&u00);
  MatDestroy(&C);
  MatDestroy(&M_sol_2);
  par.finalize();
  SlepcFinalize();
}
