/*=============================================================================
  This file is part of the code MAD
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2017/2021,
    
     Felipe Galarce

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

#include <soil.hpp>

int main(int argc, char *argv[]){

  Parameters par = MADinitialize(argc, argv);
  par.print();

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);
    
  Geometry geo;
  geo.initialize(par, io);

  Boundary boundary;
  boundary.initialize(par, geo);

  WaveComplex wave;
  wave.initialize(par, geo, boundary);

  Soil soil;
  soil.initialize(par, boundary, geo);

  Mat A = wave.assembleLHS();
  Vec b = zeros(geo.nbVertices*par.nbVariables());

  wave.setSolver();
  wave.setLHS(A);
  Vec u = zeros(par.nbVariables()*geo.nbVertices);
  
  soil.applyBC(A, b);

  MADprint("Wave: Solving linear system.\n");
  wave.solve(b, u);
    
  io.writeState(u);

  double norm_sol = norm(u);
  MADprint("Wave: Norm solution = ", norm_sol);

  //MADfinalize(wave, A, b, u, par);
}
