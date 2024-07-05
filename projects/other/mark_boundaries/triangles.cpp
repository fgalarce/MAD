/*=============================================================================
  This file is part of the code MAD
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

  ofstream file(par.dirResults() + "/triangles.scl");
  file << "scalar per element" << endl; 
  file << "part 1 "  << endl << "tria3" << endl; 
  for (int i = 0; i < geo.nbElementsBoundary(); i++){
    file << scientific << setprecision(5) << " " << (float)i;
    if ((i+1) % 6 == 0){
      file << endl;
    }
  }
  file.close();

  ofstream casefile(par.dirResults() + "/" + par.patientName() + ".case");
  casefile << "FORMAT" << endl;
  casefile << "type: ensight" << endl;
  casefile << "GEOMETRY" << endl;
  casefile << "model: 1 " + par.patientName() + ".geo" << endl;
  casefile << "VARIABLE" << endl;
  casefile << "scalar per element: 1 triangles triangles.scl" << endl;
  casefile.close();
  PetscFinalize();

}
