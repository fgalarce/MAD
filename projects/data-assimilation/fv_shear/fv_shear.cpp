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
    vector<double> data = stl(measures.measures());
    vector<double> u(nbVoxels);
    double max_u = 0.0;
    for (int idVoxel = 0; idVoxel < nbVoxels; idVoxel++){
      u[idVoxel] = sqrt( data[3*idVoxel]   * data[3*idVoxel]
                        +data[3*idVoxel+1] * data[3*idVoxel+1]
                        +data[3*idVoxel+2] * data[3*idVoxel+2]  );
      if (abs(u[idVoxel]) > max_u){
        max_u = abs(u[idVoxel]);
      }
    }
    double xi = 1.0 / max_u * PI; /* motion encoding efficiency */
    vector<double> phi_w(nbVoxels);
    vector<double> phi_u(nbVoxels);
    for (int idVoxel = 0; idVoxel < nbVoxels; idVoxel++){
      phi_u[idVoxel] = xi * u[idVoxel];
      phi_w[idVoxel] = phi_u[idVoxel] - 2*PI * floor((phi_u[idVoxel] + PI) / (2.0*PI));
    }

    /* plot */
    Vec desp = zeros(nbDofs);
    Vec phase_u = zeros(nbDofs);
    Vec phase_w = zeros(nbDofs);
    for (int i = 0; i < nbVoxels; i++){
      code = VecAXPY(desp, u[i], measures.rieszRepresenters(3*i)); CHKERR(code);
      code = VecAXPY(phase_u, phi_u[i], measures.rieszRepresenters(3*i)); CHKERR(code);
      code = VecAXPY(phase_w, phi_w[i], measures.rieszRepresenters(3*i)); CHKERR(code);
    }
    io.writeState(desp, "u", iteration);
    io.writeState(phase_u, "phi_u", iteration);
    io.writeState(phase_w, "phi_w", iteration);
  }
  measures.finalize();
  MADfinalize(par);
}
