/*=============================================================================
  This file is part of the code MAD 
  Multi-physics for mechanicAl engineering and Data assimilation
  Copyright (C) 2021-2023,
    
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

#include <brain.hpp>

int main(int argc, char *argv[]){

  Parameters par = MADinitialize(argc, argv);

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
  Vec u = zeros(io.nbVertices()*(par.nbDofsPerNode()[0] + 1));

  Mat C = ped_up.assembleLHS();
  applyBC(boundary, par, C);

  double normLHS = norm(C);
  MADprint("PEDup: norm LHS blocked = ", normLHS);

  ped_up.setSolver();
  ped_up.setLHS(C);

  ofstream file_norm_solution(par.dirResults() + "/norm_solution.txt");
  
  for (double t : range(0, par.nbIterations()*par.timeStep(), par.timeStep())){

    MADprint("\n- - - - - - - - - - - - - - - - - - - - - - - - -\n");
    MADprint("Poro: Solving poro-elasticity equation: it ", t/par.timeStep());

    /* set RHS and block with boundary conditions */
    boundary.time(t); 
    Vec b = ped_up.assembleRHS();
    applyBC(boundary, par, b);
    double normRHS = norm(b);
    MADprint("PEDup: norm RHS blocked = ", normRHS);

    ped_up.solve(b,u);

    if (t/par.timeStep() >= par.start() && ((int)floor(t/par.timeStep())) % par.jump() == 0){
      io.writeState(u, t + par.timeStep());
    }

    ped_up.update(u, t);

    double norm_sol = norm(u);
    MADprint("PEDup: Norm solution = ", norm_sol);
    file_norm_solution << norm_sol << endl; 

    VecDestroy(&b);
  }
  file_norm_solution.close();

  /* Call destructors and finish MAD */
  VecDestroy(&u);
  MatDestroy(&C);
  ped_up.finalize();
  par.finalize();
  SlepcFinalize();
}
