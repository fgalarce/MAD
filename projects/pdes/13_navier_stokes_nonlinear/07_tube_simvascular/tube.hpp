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

vector<double> inlet(vector<double> x, double t, Parameters par){
  double tstar = t/(par.distalResistances()[0]*par.capacitances()[0]);
  double area = 12.566370614;
  vector<double> bc(1, par.amplitude()/area*sin(tstar/2.0)*sin(tstar/2.0));
  return bc;
}

vector<double> noslip(vector<double> x, double t, Parameters par){
  vector<double> bc(1, 0.0);
  return bc;
}

class Tube{

  public:

    Tube(){}
    ~Tube(){};

    void initialize(Boundary & boundary, Geometry & geometry, Parameters parameters, Calculus & calculus){
      par = parameters;
      geo = geometry;
      bd = boundary;
      m_calculus = calculus;
      m_ip.initialize(par, geo, bd);
      p_d0.resize(par.resistances().size(), 0.0);
      p_d0_nonlinear.resize(par.resistances().size(), 0.0);
      p_out_nonlinear.resize(par.resistances().size(), 0.0);
      for (int i = 0; i < par.outlets().size(); i++){
        p_d0[i] = par.distalPressures0()[i];
        p_d0_nonlinear[i] = par.distalPressures0()[i];
      }

      MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
    }

    void applyBC(Vec u, Mat A, Vec b, double t){

      bd.time(t);
      bd.Dirichlet(par.inlet(), noslip, 0);
      bd.Dirichlet(par.inlet(), noslip, 1);
      bd.Dirichlet(par.inlet(), inlet, 2);
      for (int i : par.walls()){
        for (int comp = 0; comp < par.nbVariables()-1; comp++){
          bd.Dirichlet(i, noslip, comp);
        }
      }

      for (int i = 0; i < par.outlets().size(); i++){
        double Q_out = bd.flow(u, par.outlets()[i]);

        // Choose numerical integration scheme for the windkessel (0: forward Euler, 1: backward Euler, 2: RK4)
        int method = 2;
        double p_d = 0.0;
        if (method==0){
          // Semi-implicit coupling using Euler method (forward)
          p_d = p_d0[i] + par.timeStep()/par.capacitances()[i]*(Q_out - p_d0[i]/par.distalResistances()[i]);
        }
        else if (method == 1){
          // Semi-implicit coupling using Euler method (backward)
          p_d = (p_d0[i] + par.timeStep()/par.capacitances()[i]*Q_out);
          p_d /= ( 1.0 + par.timeStep()/( par.distalResistances()[i]*par.capacitances()[i] ) );
        }
        else if (method == 2){
          // Semi-implicit coupling using RK4 method with interior iterations
          double Niter = 1000;                    // number of RK iterations
          double dt_rk = par.timeStep()/Niter;  // dt for RK iterations
          double p_d0_rk = p_d0[i];             // distal pressure for RK iterations
          double k1 = 0, k2 = 0, k3 = 0, k4 = 0;
          for (int jrk=0; jrk < Niter; jrk++){
            // RK estimation
            k1 = 1/par.capacitances()[i]*(Q_out - p_d0_rk/par.distalResistances()[i]);
            k2 = 1/par.capacitances()[i]*(Q_out - (p_d0_rk + 0.5*k1*dt_rk)/par.distalResistances()[i]);
            k3 = 1/par.capacitances()[i]*(Q_out - (p_d0_rk + 0.5*k2*dt_rk)/par.distalResistances()[i]);
            k4 = 1/par.capacitances()[i]*(Q_out - (p_d0_rk + k3*dt_rk)/par.distalResistances()[i]);
            p_d = p_d0_rk + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)*dt_rk;

            // Update distal pressure
            p_d0_rk = p_d;
          }
        }
        double p_out = par.resistances()[i]*Q_out + p_d;
        p_d0_nonlinear[i] = p_d;
        p_out_nonlinear[i] = p_out;

        MADprint("  Windkessel: Boundary label: ", par.outlets()[i]);
        MADprint("  Windkessel: Flow = ", Q_out);
        MADprint("  Windkessel: Distal pressure = ", p_d);
        MADprint("  Windkessel: Proximal pressure = ", par.resistances()[i] * Q_out);
        MADprint("  Windkessel: Inwards force magnitude = ", p_out);
        bd.NeumannNormalConstant(par.outlets()[i], -p_out);

        // Add backflow stabilization
//        if (par.backflowStab()){
//          bd.backflowNeumann(u, par.outlets()[i]);
//        }
//        exit(1);

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
        Mat C = ns.assembleLHS(A, u_iter);
        b_iter = ns.assembleRHS();
        applyBC(u_iter, C, b_iter, t);

        /* solve linear system.  */
        ns.setLHS(C);

        ns.solve(b_iter, u);

        /* u - u0 */
        double rel_error = computeError(u, u_iter);
        MADprint("Navier Stokes: Non linear solver. Iteration = " + to_string(nn_iteration) + ". Residual = " + to_string(rel_error) + "\n");

        /* update field */
        VecZeroEntries(u_iter);
        VecAXPY(u_iter, 1.0, u);

        VecDestroy(&b_iter);
        if (par.nonLinearTolerance() > rel_error) break;
      }
      p_d0 = p_d0_nonlinear;
      wkssl_pressure.push_back(p_out_nonlinear[0]);
      VecDestroy(&u_iter);
      VecDestroy(&b_iter);

      // Export windkessel pressure to txt
      if (m_world_rank == 0){
        std::ofstream fout(par.dirResults() + "/wkssl_pressure.txt");
        fout << std::setprecision(10);
        for(auto const& x : wkssl_pressure)
         fout << setprecision(6) << scientific << x << '\n';
      }
    }

    double computeError(Vec u, Vec u0){
      Vec e = zeros(geo.nbVertices*par.nbVariables());
      VecAXPY(e, 1.0, u);
      VecAXPY(e, -1.0, u0);
      m_ip(e,e);
      return sqrt(m_ip(e,e)/m_ip(u,u));
    }

  private:
    Boundary bd;
    Geometry geo;
    Parameters par;
    Calculus m_calculus;
    InnerProduct m_ip;
    vector<double> p_d0;
    vector<double> p_d0_nonlinear;
    vector<double> wkssl_pressure;
    vector<double> p_out_nonlinear;
    int m_world_rank;
};
