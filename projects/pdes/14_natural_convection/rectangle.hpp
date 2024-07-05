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

vector<double> derecha(vector<double> x, double t, Parameters par){
  vector<double> bc(par.nbDofsPerNode()[2], 0.0);
  return bc;
}

vector<double> izquierda(vector<double> x, double t, Parameters par){
  vector<double> bc(par.nbDofsPerNode()[2], par.amplitude());
  return bc;
}

vector<double> noslip(vector<double> x, double t, Parameters par){
  vector<double> bc(par.nbDofsPerNode()[0], 0.0);
  return bc;
}

class FPC{

  public:

    FPC(){}
    ~FPC(){};

    void initialize(Boundary & boundary, Geometry & geometry, Parameters parameters){
      par = parameters;
      geo = geometry;
      bd = boundary;

      p_d0.resize(par.outlets().size(), 0.0);
    }

    /* update time, declare and impose boundary conditions */
    void applyBC(Mat A, Vec b, double t){
      bd.time(t); 
      bd.Dirichlet(par.inlet(), izquierda, 2);
      bd.Dirichlet(par.outlet(), derecha, 2);
      for (int i : par.walls()){
        bd.Dirichlet(i, noslip, 0);
      }
      bd.block(A);
      bd.block(b);
    }

    void solve(NaturalConvection & ns, Mat A, Boundary & boundary, Vec u, Calculus & calculus, double t, IO & io){  

      int nbDofs = ns.nbDofs;
      double TOL_nonlinear = 5e-2;
      int max_iterations = 1000;
      double rel_error = 999;
      Vec u_iter = vec(nbDofs);
      VecAXPY(u_iter, 1.0, ns.u0);
      Vec b_iter = zeros(nbDofs);
      for (int nn_iteration = 0; nn_iteration < max_iterations; nn_iteration++){
        ns.getViscosity(u_iter);
        Mat C = ns.assembleLHS(A);
        Vec b_iter = ns.assembleRHS(u_iter);
        applyBC(C, b_iter, t);

        /* solve linear system */
        ns.setLHS(C);
        ns.solve(b_iter, u);
        if (par.backflowStab()){ boundary.backflow(u); }

        /* u - u0 */
        Vec velocity = calculus.split(u)[0];
        Vec pressure = calculus.split(u)[1];
        int nbDofsVel;
        VecGetSize(velocity, &nbDofsVel);
        Vec error = zeros(nbDofsVel);
        VecAXPY(error, 1.0, velocity);
        VecAXPY(error, -1.0, calculus.split(u_iter)[0]);
        Vec error_pressure = zeros(io.nbVertices());
        VecAXPY(error_pressure, 1.0, pressure);
        VecAXPY(error_pressure, -1.0, calculus.split(u_iter)[1]);
        double rel_error = norm(error);
        double rel_error_pressure = norm(error_pressure);
        if (norm(velocity) == 0){ 
          rel_error = 0.0;
        } else {
         rel_error = rel_error/norm(velocity) + rel_error_pressure/norm(pressure);
        }
        MADprint("Navier Stokes: Non linear solver. Iteration = " + to_string(nn_iteration) + ". Residual = " + to_string(rel_error) + "\n");

        /* update field */
        VecZeroEntries(u_iter);
        VecAXPY(u_iter, 1.0, u);

        MatDestroy(&C);
        if (TOL_nonlinear > rel_error) break;

      }
    }

  private:
    Boundary bd;
    Geometry geo;
    Parameters par;
    vector<double> p_d0;
    ofstream outfile_windk;
};
