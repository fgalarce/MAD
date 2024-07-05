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

#include <mad.hpp>

int main(int argc, char *argv[]){

  Parameters par = MADinitialize(argc, argv);
  PetscErrorCode code;

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);

  Geometry geo;                 
  geo.initialize(par, io);

  InnerProduct ip; 
  ip.initialize(par, geo);

  int nbVertices = io.nbVertices();
  int nbDofs = par.nbDofsPerNode()[0]*nbVertices;
  int nbDofsPerNode = par.nbDofsPerNode()[0];

  vector<Vec> Vn_i(par.nbModes());
  vector<Vec> Vn_j(par.nbModes());
  for (int i = 0; i< par.nbModes(); i++){
    vec(Vn_i[i], nbDofs);
    vec(Vn_j[i], nbDofs);
    io.loadVector(Vn_i[i], par.space_1() + "." + wildcard(i) + ".vct");
    io.loadVector(Vn_j[i], par.space_2() + "." + wildcard(i) + ".vct");
  }

  Mat G = mat(par.nbModes(), par.nbModes());
  for (int k = 0; k < par.nbModes(); k++){
    for (int l = 0; l < par.nbModes(); l++){
      code = MatSetValue(G, k, l, ip(Vn_i[k], Vn_j[l]), INSERT_VALUES); CHKERRQ(code);
    }
  }
  MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY);
  
  vector<double> svG = getSingularValues(G, par.nbModes());

  MADprint("Distance between spaces\n    << " + par.space_1() + " >> \n and \n    " + par.space_2() + " >> \n    = ", sqrt(1.0 - min(svG)*min(svG)));
  double d2 = 1.0 - min(svG)*min(svG);
  ofstream out_file(par.dirResults() + "/value.txt");
  out_file << d2;
  out_file.close();

  MADfinalize(par);
}
