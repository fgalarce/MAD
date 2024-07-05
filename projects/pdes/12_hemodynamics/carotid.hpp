/*=============================================================================
  This file is part of the code MAD 
  Multi-physics for mechanicAl engineering and Data assimilation
  Copyright (C) 2021-2023
    
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
  if (t < 0.36){
    return -par.amplitude() * sin(PI/0.36 * t);
  } else if (t >= 0.36 and t < 0.8){
    return -PI / 0.36 * (t - 0.36) * exp(-1/70*(t - 0.36));
  } else if (t >= 0.8 and t < 1.16){
    return -par.amplitude() * sin(PI/0.36 * (t-0.8));
  } else if (t >= 1.16 and t < 1.6){
    return -PI / 0.36 * (t - 0.36 - 0.8) * exp(-1/70*(t - 0.36 - 0.8));
  } else if (t >= 1.6 and t < 1.96){
    return -par.amplitude() * sin(PI/0.36 * (t-1.6));
  } else if (t >= 1.96 and t < 2.4){
    return -PI / 0.36 * (t - 0.36 - 1.6) * exp(-1/70*(t - 0.36 - 1.6));
  } else {
    exit(1);
  }
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
      MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
      par = parameters;
      assert(par.hemodynamics() && "ERROR: Hemodynamics module should be activated to run this example.");
      geo = geometry;
      bd = boundary;
      windkFile.open(par.dirResults()+"/windkessel.txt");
      p_d0.resize(par.resistances().size());
      for (int i = 0; i < par.outlets().size(); i++){
        p_d0[i] = par.distalPressures0()[i];
        windkFile << "meanStress" << i << " p_dist" << i << " ";
      }
      windkFile << endl;
    }

    void  applyBC(Vec u, Mat A, Vec b, double t){

      bd.time(t); 
      bd.DirichletNormalParaboloid(par.inlet(), inlet);
      for (int i : par.walls()){
        bd.Dirichlet(i, noslip);
      }
      if (par.resistances().size() > 0){
        for (int i = 0; i < par.outlets().size(); i++){
          double Q_out = bd.flow(u, par.outlets()[i]);
          double p_d = 1.0/(par.capacitances()[i]/par.timeStep() + 1.0/par.distalResistances()[i]) * (Q_out + p_d0[i]*par.capacitances()[i]/par.timeStep());
          double p_out = p_d + par.resistances()[i] * Q_out;
          p_d0[i] = p_d;

          windkFile << p_out << " " << p_d << " ";
          MADprint("  Windkessel: Boundary label: ", par.outlets()[i]);
          MADprint("  Windkessel: Flow = ", Q_out);
          MADprint("  Windkessel: Distal pressure = ", p_d);
          MADprint("  Windkessel: Proximal pressure = ", par.resistances()[i] * Q_out);
          MADprint("  Windkessel: mean stress = ", p_out);
          bd.NeumannNormalConstant(par.outlets()[i], -1.0*p_out, 0);
        }
        windkFile << endl;
      }

      bd.block(A);
      bd.block(b);
    }

  private:
    ofstream windkFile;
    vector<double> p_d0;
    Boundary bd;
    Geometry geo;
    Parameters par;
    int world_rank;

};
