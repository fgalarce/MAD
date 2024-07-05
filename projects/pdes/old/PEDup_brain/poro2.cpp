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

vector<double> csf(vector<double> x, double t, Parameters par){
  vector<double> bc(1, par.csf_pressure());
  return bc;
}

//double MRE(vector<double> x, double t, Parameters par){
//  double r = 7.5;
//  double bc = par.amplitude() * ( sin(2*PI * t / par.period()) ) * ( ((r - x[1])*(r - x[1])) / (r*r) - 1.0 );
////  double bc = par.amplitude() * ( sin(2*PI * t / par.period()) );
//  if (bc > 0) {
//    bc = 0.0;
//  }
//  return bc;
//}

vector<double> MRE(vector<double> x, double t, Parameters par){
  double r = 7.5;
  vector<double> bc(2,0.0);
  bc[0] = par.amplitude() * ( sin(2*PI * t / par.period()) ) * ( 1.0 - ((r - x[1])*(r - x[1])) / (r*r) );
//  double bc = par.amplitude() * ( sin(2*PI * t / par.period()) );
  if (bc[0] > 0) {
    bc[0] = 0.0;
  }
  return bc;
}

vector<double> neck(vector<double> x, double t, Parameters par){
  vector<double> bc(3, 0.0);
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

  Mat M_sol_1 = mat(nbDofs, nbDofs); /* Mass matrix elasto */
  Mat M_sol_2 = mat(nbDofs, nbDofs); /* Mass matrix div */
  Mat Mp = mat(nbDofs, nbDofs); /* Mass pressure */
  Mat C = mat(nbDofs, nbDofs); /* Static matrix */

  assembleSystem(C, Mp, M_sol_1, M_sol_2, par, geo, fe);

//  boundary.NeumannNormal(par.bcNeumann()[0], MRE, 0);
  boundary.Dirichlet(2, MRE, 0);
  boundary.Dirichlet(par.fixed_boundary(), neck, 0);

  /* Uniform pressure at boundaries  */
  for (int i : par.walls()){
    boundary.Dirichlet(i, csf, 1);
  }

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

    /* RHS pressure */
    Vec bp = zeros(nbDofs);
    code = MatMult(Mp, u0, bp); CHKERRQ(code);

    /* RHS elastodynamics */
    code = VecScale(u0, 2.0); CHKERRQ(code);
    code = VecAXPY(u0, -1.0, u00); CHKERRQ(code);
    Vec b = zeros(nbDofs);
    code = MatMult(M_sol_1, u0, b); CHKERRQ(code);
    double Mu_sol_norm = norm(b);
    if (world_rank == 0) cout << "Poro: Norm b_solid_1 = " << Mu_sol_norm << endl;

    code = VecAXPY(b, 1.0, b_sol2); CHKERRQ(code);
    code = VecAXPY(b, 1.0, bp); CHKERRQ(code);

    /* Pulse and fixed boundary */
//    boundary.NeumannNormal(par.bcNeumann()[0], MRE, 0);
    boundary.Dirichlet(2, MRE, 0);
    boundary.Dirichlet(par.fixed_boundary(), neck, 0);
    /* Uniform pressure at boundaries  */
    for (int i : par.walls()){
      boundary.Dirichlet(i, csf, 1);
    }

   
    /* block and set up system of equations */ 
    boundary.block(b);
    if (world_rank == 0) cout << "Poro: Solving linear system" << endl;
    code = KSPSolve(ksp, b, u); CHKERRQ(code); 
    int its;
    double rnorm;
    code = KSPGetResidualNorm(ksp, &rnorm); CHKERRQ(code);
    code = KSPGetIterationNumber(ksp, &its); CHKERRQ(code);

    if (world_rank == 0) cout << "Poro: " << its << " iterations. Residual norm = " << rnorm << endl;
    
    if (t/par.timeStep() >= par.start() && ((int)floor(t/par.timeStep())) % par.jump() == 0){
      io.writeState(u, t + par.timeStep());
    }

    /* Update solutions for centered finite differences */ 
    if (t > 0){
      code = VecZeroEntries(u00); CHKERRQ(code);
      code = VecAXPY(u00, 1.0, u0); CHKERRQ(code);
    }
    code = VecZeroEntries(u0); CHKERRQ(code);
    code = VecAXPY(u0, 1.0, u); CHKERRQ(code);

//    vector<Vec> up_split = split(u, par.nbDofsPerNode(), geo.nbVertices);
    double norm_sol, norm_u, norm_p;
    norm_sol = norm(u);
//    norm_u = norm(up_split[0]);
//    norm_p = norm(up_split[1]);
//    if (world_rank == 0) cout << "Poro: Norm solution (up | u | p)  = " << norm_sol << " | " << norm_u << " | " << norm_p  << endl;
    if (world_rank == 0) cout << "Poro: Norm solution = " << norm_sol << endl;
    file_norm_solution << norm_sol << endl; 

  }
  file_norm_solution.close();

  KSPDestroy(&ksp);
  VecDestroy(&u);
  VecDestroy(&u0);
  VecDestroy(&u00);
  MatDestroy(&C);
  MatDestroy(&M_sol_1);
  MatDestroy(&M_sol_2);
  par.finalize();
  SlepcFinalize();
}
