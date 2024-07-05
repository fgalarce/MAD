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

#include<beam.hpp>

int main(int argc, char *argv[]){

  Parameters par = MADinitialize(argc, argv);

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);
    
  Geometry geo;
  geo.initialize(par, io);

  Elasticity elasticity;
  elasticity.initialize(par, geo);

  /* Boundary conditions */
  Boundary boundary;
  boundary.initialize(par, geo);

  BoundaryCondition bc;
  bc.initialize(boundary, geo, par);

  Mat A = elasticity.assembleLHS();

  int nbDofs = par.nbDofsPerNode()[0]*geo.nbVertices;
  Vec u = zeros(nbDofs);
  Vec b = zeros(nbDofs);

//  /* add source term */
//  bc.addSource(b);

  /* apply boundary conditions */ 
  bc.applyBC(A, b); 
  double normLHS = norm(A);
  double normRHS = norm(b);
  MADprint("Elasticity: Norm LHS blocked = ", normLHS);
  MADprint("Elasticity: Norm RHS blocked = ", normRHS);
  MADprint("Elasticity: Solving linear system\n");

  elasticity.setSolver();
  elasticity.setLHS(A);
  elasticity.solve(b,u);
  
  io.writeState(u, 0.0);

  double norm_sol = norm(u);
  MADprint("Elasticity: Norm solution = ", norm_sol);

  MatDestroy(&A);
  VecDestroy(&u);
  MADfinalize(par);
}
