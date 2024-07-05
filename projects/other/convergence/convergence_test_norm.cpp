/*=============================================================================
  This file is part of the code MAD
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2017/2021,
    
     Felipe Galarce at WIAS/INRIA

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

#include <mad.hpp>

int main(int argc, char *argv[]){

  /* Deploy petsc and slepc */
  PetscInitialize(&argc, &argv, NULL, NULL); 
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

  INT interpolator;
  interpolator.initialize(par, geo, io);

  /* initialize incoming geometry and field */
  IO io_i;
  io_i.initialize(par.templateGeometry(), par);

  int nbDofs = geo.nbVertices * par.nbDofsPerNode()[0];
  Vec u = vec(nbDofs);
  io.loadVector(u, par.solutionFolder1());
  Vec Iu = interpolator.interpolate_field(par.solutionFolder2(), io_i);

  Vec error = zeros(nbDofs);
  code = VecAXPY(error,  1.0,  u); CHKERRQ(code);
  code = VecAXPY(error, -1.0, Iu); CHKERRQ(code);

  double normError = norm(error);
  double normGT = norm(u);

  if (world_rank == 0) cout << "Error norm: " << normError << endl;
  if (world_rank == 0) cout << "GT norm: " << normGT << endl;

  io.writeState(Iu);

  PetscFinalize();

}
