/*=============================================================================
  This file is part of the code M.D.M.A. (or MDMA)
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

#include <ultra-4d-flow.hpp>

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

  int nbComp = geo.dimension();

  for (int i = par.start(); i < par.end(); i++){
    vector<double> field;
    vector<vector<double>> vel = importdata("/Home/flow/galarce/MAD/big_data/VEL_1_extended/vel_" + to_string(i) + ".txt");
    for (int j = 0; j < vel.size(); j++){
      for (int comp = 0; comp < nbComp; comp++){
        field.push_back(vel[j][comp]);
      }
    }
    io.writeState(field, "vel_phantom", i);
  }

  
  PetscFinalize();
}
