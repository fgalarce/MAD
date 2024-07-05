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

vector<double> pulso(vector<double> x, double t, Parameters par){
  vector<double> bc(par.nbDofsPerNode()[0], 0.0);
  bc[0] = par.amplitude();
  return bc;
}

class Soil{

  public:

    Soil(){}
    ~Soil(){};

    void initialize(Parameters parameters, Boundary & boundary, Geometry & geometry){
      par = parameters;
      geo = geometry;
      bd = boundary;
    }

    void applyBC(Mat A, Vec b){
      bd.NeumannNormalConstant(6, 100);
      bd.NeumannNormalConstant(7, 100);
      bd.Dirichlet(8, pulso, 1);
      bd.block(A);
      bd.block(b);
    }

  private:
    Boundary bd;
    Geometry geo;
    Parameters par;

};
