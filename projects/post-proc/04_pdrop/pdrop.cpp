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

  vector<vector<double>> pdrop(par.end() - par.start() + 1);

  for (int t = par.start(); t < par.end()+1; t++){
    Vec p = zeros(io.nbVertices());
    io.loadVector(p, par.dirSyntheticField() + "/" + par.variableName()[0] + "." + wildcard(t) + ".scl");
    pdrop[t-par.start()].resize(par.bdLabels().size());
    for (int i = 0; i < par.bdLabels().size(); i++){
      pdrop[t-par.start()][i] = fabs( boundary.integral(p, par.bdLabels()[i])/boundary.computeSize(par.bdLabels()[i]) - boundary.integral(p, par.inlet())/boundary.computeSize(par.inlet()) ) ;
      MADprint("Pdrop " + to_string((int)par.bdLabels()[i]) + " - " + to_string((int)par.inlet()) + ": " +  to_string(pdrop[t-par.start()][i]) + ".\n");
    }
  }

  exportData("./pdrop.txt", pdrop);
  MADfinalize(par);

}
