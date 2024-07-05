/*=============================================================================
  This file is part of the code MAD 
  Multi-physics for mechanicAl engineering and Data assimilation
  Copyright (C) 2021, 2023
    
     Felipe Galarce at INRIA/WIAS/PUCV

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

double inlet(vector<double> x, double t, Parameters par){
  return -par.amplitude()*sin(2*PI*t/par.period());
}

//vector<double> inlet(vector<double> x, double t, Parameters par){
//  vector<double> bc(3, 0.0);
//  bc[2] = par.amplitude()*sin(2*pi/0.8*t);
//  return bc;
//}

vector<double> noslip(vector<double> x, double t, Parameters par){
  vector<double> bc(3, 0.0);
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
      p_d0.resize(par.outlets().size());
      for (int i = 0; i < par.outlets().size(); i++){
        p_d0[i] = par.distalPressures0()[i];
      }
      assert(par.capacitances().size() == par.outlets().size());
      assert(par.resistances().size() == par.outlets().size());
      assert(par.distalResistances().size() == par.outlets().size());
    }

    /* update time, declare and impose boundary conditions */
    void applyBC(Vec u, Mat A, Vec b, double t){
      bd.time(t); 
      bd.DirichletNormalParaboloid(par.inlet(), inlet, 0);

      for (int i = 0; i < par.outlets().size(); i++){
        double Q_out = bd.flow(u, par.outlets()[i]);
        double p_d = 1.0/(par.capacitances()[i]/par.timeStep() + 1.0/par.distalResistances()[i]) * (Q_out + p_d0[i]*par.capacitances()[i]/par.timeStep());
        double p_out = p_d + par.resistances()[i] * Q_out;
        p_d0[i] = p_d;

        MADprint("  Windkessel: Boundary label: ", par.outlets()[i]);
        MADprint("  Windkessel: Flow = ", Q_out);
        MADprint("  Windkessel: Distal pressure = ", p_d);
        MADprint("  Windkessel: Proximal pressure = ", par.resistances()[i] * Q_out);
        MADprint("  Windkessel: Inwards force magnitude = ", p_out);
//        bd.NeumannNormalConstant(par.outlets()[i], -1.0*p_out, 0);
        bd.NeumannNormalConstant(par.outlets()[i], -1.0*p_out, 0);
      }
      for (int i : par.walls()){
        bd.Dirichlet(i, noslip, 0);
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
