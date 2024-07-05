/*=============================================================================
  This file is part of the code MAD 
  Multi-physics for mechanicAl engineering and Data assimilation
  Copyright (C) 2021-2023,
    
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

vector<double> right(vector<double> x, double t, Parameters par){
  // gradual closure -> T_star > 0
  double T_star = par.closureTime();
  vector<double> bc(par.nbDofsPerNode()[0]);
  if (t < T_star){
    bc[0] = (T_star - t)/T_star;
  } else {
    bc[0] = 0.0;
  }
  return bc;
}

vector<double> left(vector<double> x, double t, Parameters par){
  vector<double> bc(par.nbDofsPerNode()[0], 1.0);
  return bc;
}

class Pipe{

  public:

    Pipe(){}
    ~Pipe(){};

    void initialize(Parameters parameters, Boundary & boundary, Geometry & geometry){
      par = parameters;
      geo = geometry;
      bd = boundary;
    }

    /* update time, declare and impose boundary conditions */
    void applyBC(Mat A, Vec b, double t){
      bd.time(t); 
      /* Valve  */
      bd.Dirichlet(1, right, 0);
      /* Hidrostatic pressure */
      bd.Dirichlet(0, left, 1);
      bd.block(A, "symmetric");
      bd.block(b, "symmetric");
    }

    double computeError(Vec u, Vec u0){
      Vec e = zeros(geo.nbVertices*par.nbVariables());
      VecAXPY(e, 1.0, u);
      VecAXPY(e, -1.0, u0);
      double eTe, uTu;
      VecDot(e, e, &eTe);
      VecDot(u, u, &uTu);
      double error = sqrt(eTe/uTu);
      VecDestroy(&e);
      return error;
    }

    void solve(Mat A, Vec u, double t, Hammer & hammer){
      int nbDofs;
      VecGetSize(u, &nbDofs);
      Vec u_iter = zeros(nbDofs);
      VecAXPY(u_iter, 1.0, hammer.u0);
      Vec b_iter;
      for (int nn_iteration = 0; nn_iteration < par.nonLinearMaxIterations(); nn_iteration++){
        b_iter = hammer.assembleRHS(hammer.u0);
        Mat C = hammer.assembleLHS(A, u_iter);
        applyBC(C, b_iter, t);

        /* solve linear system.  */
        hammer.setLHS(C);
        hammer.solve(b_iter, u);

        /* u - u0 */
        double rel_error = computeError(u, u_iter);
        MADprint("Hammer: Non linear solver. Iteration = " + to_string(nn_iteration) + ". Residual = " + to_string(rel_error) + "\n");

        /* Update non-linear iteration */
        VecZeroEntries(u_iter);
        VecAXPY(u_iter, 1.0, u);

        if (par.nonLinearTolerance() > rel_error){
          break;
        }
      }
      VecDestroy(&u_iter);
      VecDestroy(&b_iter);
    }


  private:
    Geometry geo;
    Parameters par;
    Boundary bd;
    InnerProduct m_ip;

};
