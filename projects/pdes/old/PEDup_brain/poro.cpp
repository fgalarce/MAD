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

  Boundary boundary;
  boundary.initialize(par, geo);

  PEDup ped_up;
  ped_up.initialize(par, geo);

  /* Solution */
  Vec u = zeros(nbDofs);

//  boundary.NeumannNormal(par.bcNeumann()[0], MRE, 0);
  boundary.Dirichlet(2, MRE, 0);
  boundary.Dirichlet(par.fixed_boundary(), neck, 0);

  /* Uniform pressure at boundaries  */
  for (int i : par.walls()){
    boundary.Dirichlet(i, csf, 1);
  }

  Mat C = ped_up.assembleLHS();
  boundary.block(C);

  if (world_rank == 0) cout << "KSP: Configuring" << endl;
  KSP ksp;
  configureKSP(ksp, par.solver(), par.preconditioner(), par.monitorKSP(), par.use_solution_as_guess_KSP(), par.reuse_preconditioner());
  if (world_rank == 0) cout << "KSP: Seting operators" << endl;
  code = KSPSetOperators(ksp, C, C); CHKERR(code); 
  code = KSPSetUp(ksp); CHKERRQ(code);

  ofstream file_norm_solution(par.dirResults() + "/norm_solution.txt");
  
  for (double t : range(0, par.nbIterations()*par.timeStep(), par.timeStep())){

    if (world_rank == 0) cout << "\n- - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
    if (world_rank == 0) cout << "Poro: Solving poro-elasticity equation: it " <<  t/par.timeStep() << " t=" << t << endl;

    /* RHS Darcy */
    boundary.time(t); 

    /* Pulse and fixed boundary */
//    boundary.NeumannNormal(par.bcNeumann()[0], MRE, 0);
    boundary.Dirichlet(2, MRE, 0);
    boundary.Dirichlet(par.fixed_boundary(), neck, 0);
    /* Uniform pressure at boundaries  */
    for (int i : par.walls()){
      boundary.Dirichlet(i, csf, 1);
    }

    /* block and set up system of equations */ 
    Vec b = ped_up.assembleRHS();
    boundary.block(b);
    if (world_rank == 0) cout << "PEDup: Solving linear system" << endl;
    code = KSPSolve(ksp, b, u); CHKERRQ(code); 
    int its;
    double rnorm;
    code = KSPGetResidualNorm(ksp, &rnorm); CHKERRQ(code);
    code = KSPGetIterationNumber(ksp, &its); CHKERRQ(code);

    if (world_rank == 0) cout << "PEDup: " << its << " iterations. Residual norm = " << rnorm << endl;
    
    if (t/par.timeStep() >= par.start() && ((int)floor(t/par.timeStep())) % par.jump() == 0){
      io.writeState(u, t + par.timeStep());
    }

    /* Update solutions for centered finite differences */ 
    if (t > 0){
      code = VecZeroEntries(ped_up.u00); CHKERRQ(code);
      code = VecAXPY(ped_up.u00, 1.0, ped_up.u0); CHKERRQ(code);
    }
    code = VecZeroEntries(ped_up.u0); CHKERRQ(code);
    code = VecAXPY(ped_up.u0, 1.0, u); CHKERRQ(code);

    double norm_sol = norm(u);
    if (world_rank == 0) cout << "PEDup: Norm solution = " << norm_sol << endl;
    file_norm_solution << norm_sol << endl; 

  }
  file_norm_solution.close();

  /* Call destructors and finish MAD */
  KSPDestroy(&ksp);
  VecDestroy(&u);
  MatDestroy(&C);
  ped_up.finalize();
  par.finalize();
  SlepcFinalize();
}
