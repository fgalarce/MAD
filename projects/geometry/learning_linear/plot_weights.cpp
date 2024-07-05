/*=============================================================================
  This file is part of the code MAD
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2020,
    
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

#include<ultra-4d-flow.hpp>

int main(int argc, char *argv[]){

  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  /* Parse data file */
  assert(argc == 2);
  Parameters par(argv[1]);
  par.print();

  IO io;
  io.initialize(par); 

  for (int i = 0; i < par.nbMeshes(); i++){
    io.writeVTKgrid(importdata1D(par.dirResults() + "/weights" + wildcard(i) + ".txt"), 6, 50, 6, "weights." + wildcard(i) + ".vtk"); 
  }

  par.finalize();
  SlepcFinalize();
}
