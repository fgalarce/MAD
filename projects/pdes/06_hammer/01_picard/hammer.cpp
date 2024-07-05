/*=============================================================================
  This file is part of the code MAD
  Multi-physics for mechanicAl engineering and Data assimilation
  Copyright (C) 2017-2024,
    
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

#include <pipe.hpp>

int main(int argc, char *argv[]){

  Parameters par = MADinitialize(argc, argv);

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);
    
  Geometry geo;
  geo.initialize(par, io);

  Hammer hammer;
  hammer.initialize(par, geo);

  Boundary boundary;
  boundary.initialize(par, geo);

  Pipe pipe;
  pipe.initialize(par, boundary, geo);

  int nbDofs = (par.nbDofsPerNode()[0] + par.nbDofsPerNode()[1])*geo.nbVertices;
  Mat A = hammer.assembleLHS_static();
  Vec u = zeros(nbDofs);
  
  hammer.setSolver();

  int nbControlPoints = 3;
  vector<vector<double>> ctrlPointsU = stl_mat(par.nbIterations(), nbControlPoints);
  vector<vector<double>> ctrlPointsP = stl_mat(par.nbIterations(), nbControlPoints);
  vector<double> normSolution(par.nbIterations());
  vector<int> indexControl(nbControlPoints);
  indexControl[0] = 0;
  indexControl[1] = (int)floor((float)geo.nbVertices/2.0);
  indexControl[2] = geo.nbVertices-1;
  for (int time_iteration : range(par.nbIterations())){
    double t = time_iteration*par.timeStep();
    MADprint("\n- - - - - - - - - - - - - - - - - - - \n");
    MADprint("Solving Water-Hammer equation for time: ", t);
    MADprint("Iteration: ", time_iteration);

    vector<double> u0 = stl(hammer.u0);
    exportData(par.dirResults() + "/sol." + wildcard(time_iteration) + ".txt", u0);

    for (int i = 0; i < nbControlPoints; i++){
      ctrlPointsU[time_iteration][i] = u0[indexControl[i]];
      ctrlPointsP[time_iteration][i] = u0[indexControl[i] + io.nbVertices()];
    }
    
    /* assemble and solve */
    pipe.solve(A, u, t, hammer);

//    exportData(par.dirResults() + "/fRe." + wildcard(time_iteration) + ".txt", hammer.fRe);
//    exportData(par.dirResults() + "/Re." + wildcard(time_iteration) + ".txt", hammer.Reynolds);

    hammer.update(u, t);    

    double norm_sol = norm(u);
    normSolution[time_iteration] = norm_sol;
    MADprint("Norm solution: ", norm_sol);
  }
  exportData(par.dirResults() + "/ctrlU.txt", ctrlPointsU);
  exportData(par.dirResults() + "/ctrlP.txt", ctrlPointsP);
  exportData(par.dirResults() + "/normSolution.txt", normSolution);

  MADfinalize(hammer, par);
}
