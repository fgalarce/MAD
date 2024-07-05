/*=============================================================================
  This file is part of the code MAD 
  Multi-physics for mechanicAl engineering and Data assimilation
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

///* Functions passed as pointers for bd conditions */
//vector<double> outlet(vector<double> x, double t, Parameters par){
//  vector<double> bc(3, 0.0);
//  bc[2] = -0.1;
//  return bc;
//}

double inlet(vector<double> x, double t, Parameters par){
  return par.amplitude();
}

vector<double> noslip(vector<double> x, double t, Parameters par){
  vector<double> bc(2, 0.0);
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

    /* update time, declare and impose boundary conditions */
    void applyBC(Vec u, Mat A, Vec b, double t){
      bd.time(t); 
      bd.DirichletNormalParaboloid(par.inlet(), inlet);

      for (int i : par.walls()){
        bd.Dirichlet(i, noslip);
      }
      bd.block(A);
      bd.block(b);
    }

  private:
    Boundary bd;
    Geometry geo;
    Parameters par;
    vector<double> p_d0;

};
