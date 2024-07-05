/*=============================================================================
  This file is part of the code M.D.M.A. (or MDMA)
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2017/2021,
    
     Felipe Galarce at INRIA / WIAS

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

#include<mad.hpp>

using namespace std;

int main(int argc, char *argv[], char *envp[]){

  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  /* Parse parameters file */
  assert(argc == 2);
  Parameters par(argv[1]);
  par.print();

  IO io;
  io.initialize(par);
    
  Geometry geo;
  geo.initialize(par, io);

  CFD cfd;
  cfd.initialize(par, geo, io);

  for (int i = par.start(); i < par.end(); i++){

    if (world_rank == 0) cout << "\n- - - - - - - - - - - - - -" << endl;
    if (world_rank == 0) cout << "Iteration " << i << endl;

    /* Load velocity */
    Vec u = vec(io.nbVertices()*geo.dimension());
    loadVec(u, par.maniFolder() + "." + wildcard(i) + ".vct");

    if (world_rank == 0) cout << "Computing vorticity for u" << endl;
    Vec vorticity = cfd.vorticity(u);
    io.writeState(vorticity, "vorticity", i);

    VecDestroy(&vorticity);
    VecDestroy(&u);
  } 

  par.finalize();
  SlepcFinalize();

}
