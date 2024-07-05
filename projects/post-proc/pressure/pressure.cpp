/*=============================================================================
  This file is part of the code MAD 
  Multi-physics for mechanicAl engineering and Data assimilation
  Copyright (C) 2021,
    
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

#include<mad.hpp>

int main(int argc, char *argv[], char *envp[]){

  Parameters par = MADinitialize(argc, argv);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  PetscErrorCode code;

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);

  Geometry geo;
  geo.initialize(par, io);

  Boundary bd;
  bd.initialize(par, geo);
 
  vector<double> press(par.end()-par.start());

  for (int i = par.start(); i < par.end(); i++){

    MADprint("\n- - - - - - - - - - - - - -");
    MADprint("Iteration ", i);

    /* Load p* */
    Vec p = vec(geo.nbVertices);
    loadVec(p, par.dirSyntheticField() + wildcard(i) + ".scl");

    double pressure = bd.integral(p, 3); 
    MADprint("Pressure magnitude: ", pressure);
    press[i-par.start()] = pressure;
  } 
  exportData(press, par.dirResults() + "/press.txt");

  MADfinalize(par);

}
