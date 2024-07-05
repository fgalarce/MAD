/*=============================================================================
  This file is part of the code MAD
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2017/2022,
    
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

#include <assemble.hpp>

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

  Geometry Tgeo_i;   /* T(\Omega_i) */
  Geometry Tgeo_j;   /* T^(-1)\Omega_j */ 

  int nbVertices_i = io_i.nbVertices();
  int nbVertices_j = io_j.nbVertices();
  int nbDofs_i = par.nbDofsPerNode()[0]*nbVertices_i;
  int nbDofs_j = par.nbDofsPerNode()[0]*nbVertices_j; 
  int nbDofsPerNode = par.nbDofsPerNode()[0];

  MADprint("Piola: Harmonic extension of LDDMM template->target.\n");
  Vec HE_i = zeros(nbDofs_i);
  HE_i = harmonic_extension(par, geo_i, io_i, par.surfaceMapping(), "direct"); 
  vector<double> Tij = stl(HE_i); 
  io_i.writeState(HE_i, "HE(lddmm_ij)", 0.0);

  MADprint("Piola: Harmonic extension of LDDMM target->template.\n");
  Vec HE_j = zeros(nbDofs_j);
  HE_j = harmonic_extension(par, geo_j, io_j, par.surfaceMappingInverse(), "direct");
  vector<double> Tji = stl(HE_j); 
  io_j.writeState(HE_j, "HE(lddmm_ji)", 0.0);

  MADprint("Initializing T(Omega_i) and T^{-1}(Omega_j).");
  Tgeo_i.initialize(io_i, Tij);
  Tgeo_j.initialize(io_j, Tji);

  MADprint("Computing interpolation in geo_j of : ", par.dirSyntheticField());
  vector<Vec> v(par.nbModes());
  if (par.nbModes() == 1){
    v[0] = vec(nbDofs_i);
    io_i.loadVector(v[0], par.dirSyntheticField());
  } else {
    for(int i = 0; i < par.nbModes(); i++){
      v[i] = vec(nbDofs_i);
      io_i.loadVector(v[i], par.dirSyntheticField() + "." + wildcard(i) + ".vct");
    }
  }

  MADprint("Computing I v(\Omega_i)");
  int nbNodesPerElement = geo_j.dimension() + 1;

  MADprint("Piola: Interpolating HE in target geometry\n");
  INT interpolator;
  interpolator.initialize(par, Tgeo_j, io_j);
  Vec I_Tij = zeros(nbDofs_j);
  I_Tij = interpolator.interpolate_field(HE_i, io_i);

  /* Write interpolated volumetric extension of LDDMM mapping */
  io_j.writeState(I_Tij, "I_Tij", 0.0);

  vector<Vec> ITv = interpolator.interpolate_field(v, io_i);

//  CFD cfd;
//  cfd.initialize(par, geo_j, io_j);
  /* Piola transform */
  for (int idMode = 0; idMode < par.nbModes(); idMode++){
//    MADprint("Piola: computing Piola transform for mode #", idMode);
    io_j.writeState(ITv[idMode], "I_Tv", idMode);
//    Vec Pv = cfd.piola(ITv[idMode], I_Tij);
//    io_j.writeState(Pv, "Pv", idMode);
  }

  MADfinalize(par);
}
