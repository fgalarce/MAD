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

  Elastodynamics elastodynamics;
  elastodynamics.initialize(par, geo);

  /* Boundary conditions */
  Boundary boundary;
  boundary.initialize(par, geo);

  Soil soil;
  soil.initialize(par, boundary, geo);

  Mat A = elastodynamics.assembleLHS();
  Vec b = zeros(geo.nbVertices*par.nbDofsPerNode()[0]);

  elastodynamics.setSolver();
  elastodynamics.setLHS(A);
  Vec u = zeros(par.nbDofsPerNode()[0]*geo.nbVertices);

  for (double t : range(0, par.nbIterations()*par.timeStep(), par.timeStep())){

    MADprint("\n - - - - - - - - - - - - - - - - - - - \n");
    MADprint("Solving Elasticity equations for time: ", t);
    MADprint("Iteration: ", t/par.timeStep());

    /* block and assembly b */
    boundary.time(t); 

    b = elastodynamics.assembleRHS();
    soil.applyBC(A, b, t);

    MADprint("Elasticity: Solving linear system.\n");
    elastodynamics.solve(b,u);
    
    io.writeState(u, t + par.timeStep());
   
    elastodynamics.update(u, t);    

    double norm_sol = norm(u);
    MADprint("Elasticity: Norm solution = ", norm_sol);
  }

  MADfinalize(elastodynamics, A, b, u, par);
}
