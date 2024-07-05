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
  return -par.amplitude();
}

double outlet(vector<double> x, double t, Parameters par){
  return par.amplitude();
}

vector<double> noslip(vector<double> x, double t, Parameters par){
  vector<double> bc(1, 0.0);
  return bc;
}

vector<double> obstacle(vector<double> x, double t, Parameters par){
  vector<double> bc(1, sqrt(2)/2);
  return bc;
}

class FPC{

  public:

    FPC(){}
    ~FPC(){};

    void initialize(IO & inputOutput, Boundary & boundary, Geometry & geometry, Parameters parameters, Calculus & calculus){
      par = parameters;
      geo = geometry;
      bd = boundary;
      io = inputOutput;
      m_calculus = calculus;
      m_ip.initialize(par, geo, bd);
    }

    /* update time, declare and impose boundary conditions */
    void applyBC(Mat A, Vec b, double t){
      bd.time(t); 
      //bd.DirichletNormalParaboloid(par.inlet(), inlet, range(geo.dimension()));
      
      //Inlet, outlet
      bd.Dirichlet(1, noslip, 0);
      bd.Dirichlet(1, noslip, 1); 
      bd.Dirichlet(2, noslip, 0);
      bd.Dirichlet(2, noslip, 1);
      
      for (int i : par.walls()){
        //i-> moves betwwen wall groups (3, 4). The obstacle is 4
        for (int comp = 0; comp < par.nbVariables()-1; comp++){
          //comp-> moves between variables (ux, uy, p). We need to set ux and uy in the walls
          if (comp < 2) {
		  switch (i) {
		    case 4:
		    {
		      //Obstacle
		      bd.Dirichlet(i, obstacle, comp);
		      break;
		    }
		    default:
		    {
		      //Walls
		      bd.Dirichlet(i, noslip, comp);
		      break;
		    }
		  }
          } else {
            bd.Dirichlet(i, noslip, comp);
          }
        }
      }
      bd.block(A, "symmetric");
      bd.block(b, "symmetric");
    }

    void solve(Mat A, Boundary & bd, Vec u, double t, NS & ns){  

      Vec u_iter = zeros(ns.nbDofs);
      VecAXPY(u_iter, 1.0, ns.u0);
      Vec b_iter;

      for (int nn_iteration = 0; nn_iteration < par.nonLinearMaxIterations(); nn_iteration++){
        if (par.power_law_n() != 1.0){
          ns.getViscosity(u_iter);
        }
        if (par.writeNonLinearIterations()){
          io.writeState(u_iter, nn_iteration);
        }
        Mat C = ns.assembleLHS(A, u_iter);
        b_iter = ns.assembleRHS();
        applyBC(C, b_iter, t);

        /* solve linear system.  */
        ns.setLHS(C);

        ns.solve(b_iter, u);
        if (par.backflowStab()){ bd.backflow(u); }

        /* u - u0 */
        double rel_error = computeError(u, u_iter);
        MADprint("Navier Stokes: Non linear solver. Iteration = " + to_string(nn_iteration) + ". Residual = " + to_string(rel_error) + "\n");

        /* update field */
        VecZeroEntries(u_iter);
        VecAXPY(u_iter, 1.0, u);

        VecDestroy(&b_iter);
        if (par.nonLinearTolerance() > rel_error) break;
      }
      VecDestroy(&u_iter);
    }

    double computeError(Vec u, Vec u0){
      Vec e = zeros(geo.nbVertices*par.nbVariables());
      VecAXPY(e, 1.0, u);
      VecAXPY(e, -1.0, u0);
      m_ip(e,e);
      double error = sqrt(m_ip(e,e)/m_ip(u,u));
      VecDestroy(&e);
      return error;
    }

  private:
    Boundary bd;
    Geometry geo;
    Parameters par;
    Calculus m_calculus;
    IO io;
    InnerProduct m_ip;
};
