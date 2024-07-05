/*=============================================================================
  This file is part of the code M.D.M.A. (or MDMA)
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2021,
    
     Felipe Galarce at INRIA

  MDMA is free software; you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by the Free
  Software Foundation; either version 2.1 of the License, or (at your option)
  any later version.

  MDMA is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with MDMA. If not, see http://www.gnu.org/licenses/.
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

  int bdLabel = 3;
  ofstream p_file(par.dirResults() + "/p_gt.txt");
  ofstream p_file_star(par.dirResults() + "/p_star.txt");

  Vec uno = vec(geo.nbVertices);
  VecSet(uno, 1.0);
  double surface_area = boundary.integral(uno, bdLabel);
  cout << "  Surface area : " << surface_area << endl;

  for (int i = par.start(); i < par.start() + par.nbSimulations(); i++){
    MADprint("-- sim " + wildcard(i) + " --\n");
    double int_bd_avg = 0.0;
    double int_bd_avg_star = 0.0;
    for (int snapId = 0; snapId < par.snaps_per_sim(); snapId++){
      Vec p = zeros(geo.nbVertices);
      Vec q = zeros(geo.nbVertices);
//      cout << norm(p) << endl;
      loadVec(p, par.maniFolder() + "/sim" + wildcard(i) + "/pressure_GT." + wildcard(snapId) + ".scl");
      loadVec(q, par.maniFolder() + "/sim" + wildcard(i) + "/pressure_star." + wildcard(snapId) + ".scl");
      double int_bd = boundary.integral(p, bdLabel); 
      double int_bd_star = boundary.integral(q, bdLabel); 
      int_bd_avg = int_bd_avg + int_bd;
      int_bd_avg_star = int_bd_avg_star + int_bd_star;
    }
    p_file << int_bd_avg / par.snaps_per_sim() << endl;
    p_file_star << int_bd_avg_star / par.snaps_per_sim() << endl;
    
  }
  p_file.close();
  MADfinalize(par);
}
