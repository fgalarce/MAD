/*=============================================================================
  This file is part of the code MAD
  Multy-physics and Data assimilation for Medical Applications.
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

  Parameters par = MADinitialize(argc, argv);

  /* Initialize MAD objects */
  IO io_i;
  io_i.initialize(par.templateGeometry(), par.dirResults() + "/template/", par);

  IO io_j;
  io_j.initialize(par);

  /* \Omega_i */
  Geometry geo_i;
  geo_i.initialize(io_i);     

  /* \Omega_j */
  Geometry geo_j;                 
  geo_j.initialize(par, io_j);

  int nbVertices_i = io_i.nbVertices();
  int nbVertices_j = io_j.nbVertices();
  int nbDofs_i = par.nbDofsPerNode()[0]*nbVertices_i;
  int nbDofs_j = par.nbDofsPerNode()[0]*nbVertices_j; 
  int nbDofsPerNode = par.nbDofsPerNode()[0];

  MADprint("Computing I v(\Omega_i)");
  Vec v = vec(nbDofs_i);
  io_i.loadVector(v, par.dirSyntheticField());

  INT interpolator;
  interpolator.initialize(par, geo_j, io_j);
  Vec ITv = interpolator.interpolate_field(v, io_i);
  double norm_itv;
  VecNorm(ITv, NORM_2, &norm_itv);
  MADprint("norm ITv ", norm_itv);
  io_j.writeState(ITv, "ITv");

  MADfinalize(par);
}
