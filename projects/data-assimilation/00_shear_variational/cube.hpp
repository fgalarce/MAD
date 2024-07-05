/*=============================================================================
  This file is part of the code MAD 
  Multi-physics for mechanicAl engineering and Data assimilation
  Copyright (C) 2023,
    
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

Vec computeSource(Parameters par, Geometry & geo, double t){
  Vec source = zeros(geo.nbVertices*par.nbDofsPerNode()[0]);
  for (int i = 0; i < geo.nbVertices; i++){
    double xx = geo.coordinates()[i][0] - 0.5;
    double yy = geo.coordinates()[i][1] - 0.5;
    double zz = geo.coordinates()[i][2] - 1.5;
    double timeComponent = sin(2*PI/par.period()*t);
    if(xx*xx + yy*yy + zz*zz < 0.1){
      for (int j = 0; j < par.nbDofsPerNode()[0]; j++){
        vecSet(source, par.nbDofsPerNode()[0]*i + j, par.amplitude()*timeComponent);
      }
    }
  }
  VecAssemblyBegin(source);
  VecAssemblyEnd(source);
  return source;
}
