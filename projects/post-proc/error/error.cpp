/*=============================================================================
  This file is part of the code M.D.M.A. (or MDMA)
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2020,
    
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

#include <ultra-4d-flow.hpp>

int main(int argc, char *argv[]){

  assert(argc == 2);

  string data_file = argv[1];
  
  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  Parameters par(data_file);
  par.print();

  /* Initialize MDMA objects */
  IO io;
  io.initialize(par);

  Geometry geo;
  geo.initialize(par, io);

  Vec gt = vec(geo.nbVertices*3);
  loadVec(gt, par.maniFolder() + "/sim00000/wss_gt.00010.vct");

  io.writeState(gt, "gt");

  Vec app = vec(geo.nbVertices*3);
  loadVec(app, par.maniFolder() + "/sim00000/wss.00010.vct");

  io.writeState(app, "app");
  VecAXPY(gt, -1.0, app);

  io.writeState(gt, "error");

  SlepcFinalize();
}
