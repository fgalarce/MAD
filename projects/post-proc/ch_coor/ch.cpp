/*=============================================================================
  This file is part of the code MAD
  Copyright (C) 2021-2023,
    
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

  vector<double> x = range(0,1+0.1,0.1);
  vector<double> y = range(0,1+0.1,0.1);
  vector<double> z = range(0,0.5+0.1,0.1);
  vector<double> data(nx*ny*nz);
  int nx = x.size();
  int ny = y.size();
  int nz = z.size();

  io.writeVTKgrid(data, nx, ny, nz, "vtkFile");

    
  MADfinalize(par);
}
