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

  int nbDofs = io.nbVertices();
  Vec indicatrix = zeros(nbDofs);
  double * values = new double[nbDofs];
  for (int i = 0; i < nbDofs; i++){
    values[i] = (double)geo.coordinatesLabels()[i];
  }
  code = VecSetValues(indicatrix, nbDofs, &range(nbDofs)[0], values, INSERT_VALUES); CHKERRQ(code);
  code = VecAssemblyBegin(indicatrix); CHKERRQ(code);
  code = VecAssemblyEnd(indicatrix); CHKERRQ(code);

  io.writeState(indicatrix, "indicatrix");

  PetscFinalize();

}
