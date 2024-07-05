/*=============================================================================
  This file is part of the code MAD
  Multy-physics and dAta assimilation for meDical applications.
  Copyright (C) 2017-2022,
    
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

#include <ultra-4d-flow.hpp>

int main(int argc, char *argv[]){

  assert(argc == 2);

  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  Parameters par(argv[1]);
  par.print();

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);

  Geometry geo;
  geo.initialize(par, io);

  Learning learning;
  learning.initialize(par, geo);

  learning.voxelize();
  
  learning.write_geo_unstructured_grid();

  exportData(par.dirResults() + "/g.txt", learning.voxelized_geometry());

  SlepcFinalize();

}
