/*=============================================================================
  This file is part of the code MAD
  Multi-physics for mechanicAl engineering and Data assimilation
  Copyright (C) 2017-2023,
    
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

#include <toy.hpp>

int main(int argc, char *argv[]){

  Parameters par = MADinitialize(argc, argv);

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);
    
  Geometry geo;
  geo.initialize(par, io);

  Boundary boundary;
  boundary.initialize(par, geo);

  Heat heat;
  heat.initialize(par, geo, boundary);

  BoundaryCondition bc;
  bc.initialize(boundary, geo, par);

  int nbDofs = par.nbDofsPerNode()[0]*geo.nbVertices;
  Mat A = heat.assembleLHS();
  bc.applyBC(A);
  heat.setSolver();
  heat.setLHS(A);
  Vec u = zeros(nbDofs);
  Vec b = zeros(nbDofs);

  for (double t : range(0, par.nbIterations()*par.timeStep(), par.timeStep())){

    int iteration = (int)round(t/par.timeStep());
    MADprint("\n- - - - - - - - - - - - - - - - - - - \n");
    MADprint("Solving Heat equation for time: ", t);
    MADprint("Iteration: ", iteration);
    b = heat.assembleRHS();
    bc.applyBC(b, t);

    heat.solve(b,u);

    if (iteration >= par.start()){
      io.writeState(u, t + par.timeStep());
    }

    heat.update(u, t);    

    double norm_sol = norm(u);
    MADprint("Heat: Norm solution = ", norm_sol);
  }

  MADfinalize(heat, par, A, b, u);
}
