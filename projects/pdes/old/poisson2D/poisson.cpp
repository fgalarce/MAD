/*=============================================================================
  This file is part of the code AsimovHomework
  Copyright (C) 2020,
    
     Felipe Galarce at INRIA

  AsimovHomework is free software; you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by the Free
  Software Foundation; either version 2.1 of the License, or (at your option)
  any later version.

  AsimovHomework is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with AsimovHomework. If not, see http://www.gnu.org/licenses/.
  =============================================================================*/

#include <asimov_homework.hpp>
#include <poisson.hpp>

int main(int argc, char *argv[]){

  assert(argc == 2);
  
  /* Deploy petsc and slepc */
  PetscInitialize(&argc, &argv, NULL, NULL); 
 
  Parameters par(argv[1]);
  par.print();

  /* Initialize AsimovHomework */
  IO io;
  io.square_mesh(par);
    
  Geometry geo;
  geo.initialize(io);

  Poisson poisson;
  poisson.initialize(par, geo);

  /* Assemble system  */
  poisson.assembleStiffness();
  poisson.assembleSource();

  /* Impose 0.0 boundary conditions */
  poisson.applyBC(0.0);

  /* Set up Krylov solver and run GMRES+ASM  */
  Vec solution = poisson.solve();

  /* Write solution to ensight */
  io.writeState(solution);

  par.finalize();
  PetscFinalize();
}
