/*=============================================================================
  This file is part of the code MAD
  Copyright (C) 2021-2024,
    
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
  along with MDMA. If not, see http://www.gnu.org/licenses/.
  =============================================================================*/

#include <square.hpp>

int main(int argc, char *argv[]){

  MADinitialize(argc, argv);

  /* Parse data file */
  string data_file = argv[1]; 
  Parameters par(data_file);

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);
    
  Geometry geo;
  geo.initialize(par, io);

  Boundary boundary;
  boundary.initialize(par, geo);

  BoundaryCondition bc;
  bc.initialize(boundary, geo, par);

  Laplacian laplacian;
  laplacian.initialize(par, geo, boundary);

  int nbDofs = par.nbVariables()*io.nbVertices();

  Mat K = laplacian.assembleLHS(); 

  /* Assemble source */
  Vec rhs = vec(nbDofs);
  vector<double> source(nbDofs);
  for(int i = 0; i < source.size(); i++){
    source[i] = geo.coordinates()[i][0];
  }
  MatMult(laplacian.mass(), petsc(source), rhs);

  /* impose boundary conditions */
  bc.applyBC(K, rhs);

  if (nbDofs < 101){
    cout << geo.coordinates() << endl;
    MatView(K, PETSC_VIEWER_STDOUT_WORLD);
    VecView(rhs, PETSC_VIEWER_STDOUT_WORLD);
  }

  /* solver */
  Vec u = zeros(nbDofs);
  laplacian.setSolver();
  laplacian.setLHS(K);
  laplacian.solve(rhs, u);
  
  io.writeState(u);

  laplacian.finalize();
  MADfinalize(par);
}
