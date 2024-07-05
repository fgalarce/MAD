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

vector<double> empotrado(vector<double> x, double t, Parameters par){
  vector<double> bc(par.nbDofsPerNode()[0], 0.0);
  return bc;
}

class BoundaryCondition{

  public:

    BoundaryCondition(){}
    ~BoundaryCondition(){};

    void initialize(Boundary & boundary, Geometry & geometry, Parameters parameters){
      par = parameters;
      geo = geometry;
      bd = boundary;
    }

    void applyBC(Mat A, Vec b){
      for (int i=0; i<par.walls().size(); i++){
        bd.Dirichlet(par.walls()[i], empotrado);
      }
      bd.block(A);
      bd.block(b);
    }

    void addSource(Vec b){
      Vec source = zeros(geo.nbVertices*par.nbDofsPerNode()[0]);
      for (int i = 0; i < geo.nbVertices; i++){
        double source = 0;
        double xx = geo.coordinates()[i][0] - 0.5;
        double yy = geo.coordinates()[i][1] - 0.5;
        double zz = geo.coordinates()[i][2] - 1.5;
        if(xx*xx + yy*yy + zz*zz < 0.1){
          for (int j = 0; j < par.nbDofsPerNode()[0]; j++){
            vecSet(b, par.nbDofsPerNode()[0]*i + j, par.amplitude());
          }
        }
      }
      VecAssemblyBegin(b);
      VecAssemblyEnd(b);
    }

  private:
    Boundary bd;
    Geometry geo;
    Parameters par;

};
