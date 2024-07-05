/*=============================================================================
  This file is part of the code MAD 
  Multi-physics for mechanicAl engineering and Data assimilation
  Copyright (C) 2024,
    
     Felipe Galarce

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

vector<double> zero(vector<double> x, double t, Parameters par){
  vector<double> bc(1, 0.0);
  return bc;
}
vector<double> csf(vector<double> x, double t, Parameters par){
  vector<double> bc(1, par.csf_pressure());
  return bc;
}

vector<double> empotrado(vector<double> x, double t, Parameters par){
  vector<double> bc(par.nbDofsPerNode()[0], 0.0);
  return bc;
}

vector<double> MRE(vector<double> x, double t, Parameters par){
  vector<double> bc(3, 0.0); 
  bc[0] = -par.amplitude();
  return bc;
}

class Cube{

  public:

    Cube(){}
    ~Cube(){};

    void initialize(Boundary & boundary, Parameters parameters){
      par = parameters;
      bd = boundary;
    }

    void applyBC(Mat A, Vec b){

      /* = = = = Bottom = = = = */
      /* condition on pressure */
      bd.Dirichlet(par.bottom(), csf, 2);
      bd.Dirichlet(par.bottom(), zero, 3);

      /* Neumann boundary condition */
      bd.Neumann(par.bottom(), MRE, 0);

      /* = = = = Top = = = = */
      /* condition on displacements */
      bd.Dirichlet(par.top(), empotrado, 0);
      bd.Dirichlet(par.top(), empotrado, 1);

      bd.block(A, "symmetric");
      bd.block(b, "symmetric");
    }

  private:
    Boundary bd;
    Geometry geo;
    Parameters par;

};
