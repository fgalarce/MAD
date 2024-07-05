/*=============================================================================
  This file is part of the code MAD
  Copyright (C) 2017 - 2021,
    
     Felipe Galarce at WIAS

  MAD is free software; you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by the Free
  Software Foundation; either version 2.1 of the License, or (at your option)
  any later version.

  MAD is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with MDMA. If not, see http://www.gnu.org/licenses/.
  =============================================================================*/

#include <mad.hpp>

int main(int argc, char *argv[]){

  Parameters par = MADinitialize(argc, argv);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  PetscErrorCode code;

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);

  Geometry geo;
  geo.initialize(par, io);

  InnerProduct ip;
  ip.initialize(par, geo);

  LinearAlgebra la;
  la.initialize(par, ip);

  Measures measures;
  measures.initialize(par, geo, la, io);

  int m = measures.nbMeasures();
  int nbVoxels = measures.nbVoxels();

  int nbVertices = geo.nbVertices;
  int nbDofs = nbVertices*3; 
  
  for (int iteration = par.start(); iteration < par.end(); iteration++ ){

    if (world_rank == 0) cout << "\n - - - - - - - - - - - - - - - - - - - - - " << endl;
    if (world_rank == 0) cout << " Computing phase: " << iteration<< endl;

    measures.computeSynthetic(iteration);

//    Vec Measure = zeros(nbDofs);
//    for (int i = 0; i < m; i++){
//      VecAXPY(Measure, stl(measures.measures())[i], measures.rieszRepresenters()[i]);
//    }
//    io.writeState(Measure, "measures", iteration);
//    vector<vector<vector<double>>> image2d = m
//    if (world_rank == 0 ){
//      exportData(par.dirResults() + "/measures_x." + wildcard(iteration) + ".txt", image2d[0]);
//      exportData(par.dirResults() + "/measures_y." + wildcard(iteration) + ".txt", image2d[1]);
//      exportData(par.dirResults() + "/measures_z." + wildcard(iteration) + ".txt", image2d[2]);
//    }
    exportData(par.dirResults() + "/measures." + wildcard(iteration) + ".txt", stl(measures.measures()));
  }
  measures.finalize();
  MADfinalize(par);
}
