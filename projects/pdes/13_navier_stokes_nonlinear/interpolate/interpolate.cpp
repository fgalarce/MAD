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

  IO io_gt;
  io_gt.initialize(par.geometry(), "./results_int/");

  Geometry geo_gt;
  geo_gt.initialize(io_gt);

  INT interpolator;
  interpolator.initialize(par, geo_gt, io_gt);

  for (int t : range(par.start(), par.end())){

    MADprint("\n- - - - - - - - - - - - - - - - - - - - - - - - -\n");
    MADprint("Interpolating : ", t);

    Vec u = vec(par.nbDofsPerNode()[0]*io.nbVertices());
    io.loadVector(u, par.dirSyntheticField() + "/" + par.variableName()[0] + "." + wildcard(t) + ".vct");
    io.writeState(u, "velocity", t);
    Vec uInt = interpolator.interpolate_field(u, io, 0); 
    io_gt.writeState(uInt, "velocity", t);

  }

  MADfinalize(par);
}
