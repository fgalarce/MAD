/*=============================================================================
  This file is part of the code MAD 
  Multi-physics for mechanicAl engineering and Data assimilation
  Copyright (C) 2021,
    
     Felipe Galarce at INRIA

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
  bc[1] = par.amplitude() * sin(2.0*PI/par.period() * t) * (x[1] + 10.5211)/(16.86762);
  return bc;
}

vector<double> neck(vector<double> x, double t, Parameters par){
  vector<double> bc(3, 0.0);
  return bc;
}

void applyBC(Boundary boundary, Parameters par, Mat A){
  for (int i : par.walls()){
    boundary.Dirichlet(i, csf, 1);
  }
  boundary.Dirichlet(3, csf_ventricles, 1);
  boundary.Dirichlet(par.fixed_boundary(), neck, 0);
  boundary.block(A);
}

void applyBC(Boundary boundary, Parameters par, Vec b){
  boundary.Neumann(par.bcNeumann()[0], MRE, 0);
  boundary.Neumann(par.bcNeumann()[1], MRE, 0);
  for (int i : par.walls()){
    boundary.Dirichlet(i, csf, 1);
  }
  boundary.Dirichlet(3, csf_ventricles, 1);
  boundary.Dirichlet(par.fixed_boundary(), neck, 0);
  boundary.block(b);
}
