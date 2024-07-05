/*=============================================================================
  This file is part of the code MAD
  Copyright (C) 2017-2023,
    
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

#include <cube.hpp>

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

  PoroComplex poro_complex; 
  poro_complex.initialize(par, geo, boundary, io);

  poro_complex.setSolver();
  int nbDofs = 0;
  for (int i = 0; i < par.nbDofsPerNode().size(); i++){
    nbDofs += io.nbVertices()*par.nbDofsPerNode()[i];
  }

  Vec u = zeros(nbDofs);

  Cube solver;
  solver.initialize(boundary, par);

  Mat A = poro_complex.assembleLHS(); 

  Vec b = poro_complex.assembleRHS();
  solver.applyBC(A, b);

  poro_complex.setLHS(A);
  poro_complex.solve(b, u);

  double normSol = norm(u);
  MADprint("Norm l2 solution: ", normSol);

  io.writeState(u, 0);


  poro_complex.finalize();
  MADfinalize(par);
}
