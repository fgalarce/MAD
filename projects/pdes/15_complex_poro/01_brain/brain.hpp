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

vector<double> csf(vector<double> x, double t, Parameters par){
  vector<double> bc(1, par.csf_pressure());
  return bc;
}

vector<double> csf_ventricles(vector<double> x, double t, Parameters par){
  vector<double> bc(1, par.csf_pressure()*par.csf_factor());
  return bc;
}

vector<double> MRE(vector<double> x, double t, Parameters par){
  vector<double> bc(3, 0.0); 
  bc[1] = par.amplitude();
  return bc;
}

vector<double> neck(vector<double> x, double t, Parameters par){
  vector<double> bc(3, 0.0);
  return bc;
}

vector<double> zero(vector<double> x, double t, Parameters par){
  vector<double> bc(3, 0.0);
  return bc;
}

class Brain{

  public:

    Brain(){}
    ~Brain(){};

    void initialize(Boundary & boundary, Parameters parameters){
      par = parameters;
      bd = boundary;
    }

    /* update time, declare and impose boundary conditions */
    void applyBC(Mat A, Vec b){
      for (int i : par.walls()){
        bd.Dirichlet(i, csf, 2);
        bd.Dirichlet(i, zero, 3);
      }
      bd.Dirichlet(3, csf_ventricles, 2);
      bd.Dirichlet(3, zero, 3);
      bd.Dirichlet(par.fixed_boundary(), neck, 0);
      bd.Dirichlet(par.fixed_boundary(), neck, 1);
      bd.Neumann(par.bcNeumann()[0], MRE, 0);
      bd.block(A, "symmetric");
      bd.block(b, "symmetric");
    }
    
  private:
    Boundary bd;
    Parameters par;
    PetscErrorCode code;
};
