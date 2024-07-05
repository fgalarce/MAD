/*=============================================================================
  This file is part of the code MAD
  Copyright (C) 2021-2023,
    
     Felipe Galarce at INRIA/WIAS/PUCV 

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

  Parameters par = MADinitialize(argc, argv);

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);
    
  Geometry geo;
  geo.initialize(par, io);

  Boundary boundary;
  boundary.initialize(par, geo);

  Calculus calculus;
  calculus.initialize(par, geo, boundary);

  int nbDofs = io.nbVertices()*par.nbVariables();

  vector<Vec> u(geo.dimension());
  vec(u[0], io.nbVertices());
  vec(u[1], io.nbVertices());
  vec(u[2], io.nbVertices());

  for (int t : range(par.start(), par.end(), par.jump())){

    MADprint("\n- - - - - - - - - - - - - - - - - - - - - - - - -\n");
    MADprint(" Computing WSS: ", t);
    io.loadVector(u[0], par.dirSyntheticField() + "/ux." + wildcard(t) + ".scl");
    io.loadVector(u[1], par.dirSyntheticField() + "/uy." + wildcard(t) + ".scl");
    io.loadVector(u[2], par.dirSyntheticField() + "/uz." + wildcard(t) + ".scl");
    Vec wss = calculus.wallShearStress(u,par.walls()[0]);
    io.writeState(wss, "wss", t);

  }
  MADfinalize(par);
}
