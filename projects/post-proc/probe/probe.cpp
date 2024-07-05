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

#include <mad.hpp>

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

  int nbDofsPress = io.nbVertices();
  int nbVertices = io.nbVertices();

  /* locate probe nearest neighbourd */
  vector<double> proble_location(geo.dimension());
  proble_location[0] = 0.5;
  proble_location[1] = 0.15;

  int nn = -1;
  double distance = 99999;
  for (int i = 0; i < nbVertices; i++){
    double candidate_distance = norm(geo.coordinates()[i] - proble_location);
    if (candidate_distance < distance ){
      distance = candidate_distance;
      nn = i;
    }
  }

  /* Get field at probe location  */
  vector<vector<double>> velocity(par.nbIterations());
  vector<double> pressure(par.nbIterations());
  for (int i = 0; i < par.nbIterations(); i++){
    Vec vel = zeros(nbVertices*3);
    Vec pre = zeros(nbDofsPress);
    loadVec(vel, par.maniFolder() + "/" + par.variableName()[0] + "." + wildcard(i) + ".vct");
    loadVec(pre, par.maniFolder() + "/" + par.variableName()[1] + "." + wildcard(i) + ".scl");
    velocity[i].resize(3);
    for (int comp = 0; comp < 3; comp++){ /* ensight compels three vector components even in 2d */
      int index = 3*nn+comp;
      code = VecGetValues(vel, 1, &index, &velocity[i][comp]); CHKERRQ(code);
    }
    code = VecGetValues(pre, 1, &nn, &pressure[i]); CHKERRQ(code);
  }
  exportData(par.dirResults() + "/velocity_at_probe.txt", velocity);
  exportData(par.dirResults() + "/pressure_at_probe.txt", pressure);

  par.finalize();
  SlepcFinalize();
}
