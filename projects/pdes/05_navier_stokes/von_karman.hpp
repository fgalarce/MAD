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

/* Functions passed as pointers for bd conditions */
vector<double> inlet(vector<double> x, double t, Parameters par){
  vector<double> center(2, 0.0);
  vector<double> bc(2, 0.0);
  double radius = 0.25;
  center[1] = radius;
  bc[0] = par.amplitude() * (1.0 - dot(center - x, center - x) / (radius*radius));
  return bc;
}

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

      p_d0.resize(par.outlets().size(), 0.0);
      assert(par.capacitances().size() == par.outlets().size());
      assert(par.resistances().size() == par.outlets().size());
      assert(par.distalResistances().size() == par.outlets().size());

    }

    /* update time, declare and impose boundary conditions */
    void applyBC(Vec u, Mat A, Vec b, double t){
      bd.time(t); 
      bd.Dirichlet(par.inlet(), inlet);
      for (int i : par.walls()){
        bd.Dirichlet(i, noslip);
      }
//      for (int i = 0; i < par.outlets().size(); i++){
//        double Q_out = bd.flow(u, par.outlets()[i]);
//        double p_d = 1.0/(par.capacitances()[i]/par.timeStep() + 1.0/par.distalResistances()[i]) * (Q_out + p_d0[i]*par.capacitances()[i]/par.timeStep());
//        double p_out = p_d + par.resistances()[i] * Q_out;
//        p_d0[i] = p_d;
//        outfile_windk.open(par.dirResults() + "./windkessel.txt", ios_base::app); 
//        outfile_windk << Q_out << " " << p_d << " " << p_out << endl;
//
//        MADprint("  Windkessel: Boundary label: ", par.outlets()[i]);
//        MADprint("  Windkessel: Flow = ", Q_out);
//        MADprint("  Windkessel: Distal pressure = ", p_d);
//        MADprint("  Windkessel: Proximal pressure = ", par.resistances()[i] * Q_out);
//        MADprint("  Windkessel: Inwards force magnitude = ", p_out);
//        bd.NeumannNormalConstant(par.outlets()[i], -1.0*p_out, 0);
//      }
      bd.block(A, "symmetric");
      bd.block(b, "symmetric");
    }

  private:
    Boundary bd;
    Geometry geo;
    Parameters par;
    vector<double> p_d0;
    ofstream outfile_windk;
};
