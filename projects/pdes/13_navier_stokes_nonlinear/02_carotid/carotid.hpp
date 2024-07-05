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
  int iteration = (int)round(t/par.timeStep());
  double umax = par.amplitude() * importdata1D(par.inflow_data(), "shut_up")[iteration]; // 0.78367 is the inlet size
  return -umax;
}

vector<double> noslip(vector<double> x, double t, Parameters par){
  vector<double> bc(par.nbDofsPerNode()[0], 0.0);
  return bc;
}

class Aorta{

  public:

    Aorta(){}
    ~Aorta(){};

    void initialize(NS & navier_stokes, Parameters parameters, Boundary & boundary, Geometry & geometry, Calculus & calculus, IO & io){
      MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
      par = parameters;
      assert(par.hemodynamics() && "ERROR: Hemodynamics module should be activated to run this example.");
      geo = geometry;
      bd = boundary;
      m_calculus = calculus;
      m_io = io;
      ns = navier_stokes;
      p_d0.resize(par.resistances().size(), 0.0);
      p_d0_nonlinear.resize(par.resistances().size(), 0.0);
      p_d_nonlinear.resize(par.resistances().size(), 0.0);
      wkssl.resize(par.resistances().size(), 0);
      for (int i = 0; i < par.outlets().size(); i++){
        p_d0[i] = par.distalPressures0()[i];
        p_d0_nonlinear[i] = par.distalPressures0()[i] + 1e-10;
        p_d_nonlinear[i] = par.distalPressures0()[i];
      }
      m_ip.initialize(par, geo, bd);

      /* Checking compatiblity with inflow data */
      double dt_file = importdata1D(par.inflow_data(), "shut_up")[0];
      if (dt_file != par.timeStep()){
        MADprint("ERROR: inflow data file not compatible with parameters. dt_file = " + to_string(dt_file) + ". dt_parameters = " + to_string(par.timeStep()) + ".");
        exit(1);
      }

      MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
    }

    double computeError(Vec u, Vec u0){
      Vec e = zeros(geo.nbVertices*par.nbVariables());
      VecAXPY(e, 1.0, u);
      VecAXPY(e, -1.0, u0);
      double ip_error = m_ip(e,e);
//      if (ip_error == 0){
//        return 0;
//      }
      double error = sqrt(ip_error/m_ip(u,u));
      VecDestroy(&e);
      return error;
    }

    void applyBC(Vec u, Mat A, Vec b, double t){

      bd.time(t);
      bd.DirichletNormalParaboloid(par.inlet(), inlet, range(geo.dimension()));
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
          // coupling using Euler method (forward)
          p_d = p_d0[i] + par.timeStep()/par.capacitances()[i]*(Q_out - p_d0[i]/par.distalResistances()[i]);
        }
        else if (method == 1){
          // coupling using Euler method (backward)
          p_d = (p_d0[i] + par.timeStep()/par.capacitances()[i]*Q_out);
          p_d /= ( 1.0 + par.timeStep()/( par.distalResistances()[i]*par.capacitances()[i] ) );
        }
        else if (method == 2){
          // coupling using RK4 method with interior iterations
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
        p_d_nonlinear[i] = p_d;
        wkssl[i] = p_out;

        MADprint("  Windkessel: Boundary label: ", par.outlets()[i]);
        MADprint("  Windkessel: Flow = ", Q_out);
        MADprint("  Windkessel: Distal pressure = ", p_d);
        MADprint("  Windkessel: Proximal pressure = ", par.resistances()[i] * Q_out);
        MADprint("  Windkessel: Inwards force magnitude = ", p_out);
        bd.NeumannNormalConstant(par.outlets()[i], -p_out);

//        // Add backflow stabilization
//        if (par.backflowStab()){
//          bd.backflowNeumann(u, par.outlets()[i]);
//        }

      }

      bd.block(A, "symmetric");
      bd.block(b, "symmetric");
    }

    void solve(Mat A, Boundary & bd, Vec u, double t, NS & ns){

      int nbDofs = ns.nbDofs;
      Vec u_iter = zeros(nbDofs);
      VecAXPY(u_iter, 1.0, ns.u0);
      Vec b_iter;
      for (int nn_iteration = 0; nn_iteration < par.nonLinearMaxIterations(); nn_iteration++){
        if (par.power_law_n() != 1.0){
          ns.getViscosity(u_iter);
        }
        if (par.timeIntegration() == "BDF2"){
          b_iter = ns.assembleRHS(ns.u0, ns.u00);
        } else {
          b_iter = ns.assembleRHS(ns.u0);
        }
        Mat C = ns.assembleLHS(A, u_iter);
        applyBC(u_iter, C, b_iter, t);

        /* solve linear system.  */
        ns.setLHS(C);
        ns.solve(b_iter, u);

        /* u - u0 */
        double rel_error = computeError(u, u_iter);
        MADprint("Navier Stokes: Non linear solver. Iteration = " + to_string(nn_iteration) + ". Residual = " + to_string(rel_error) + "\n");

        // Check windkessels convergence
        double rel_werror = 0.0;
        for (int k = 0; k < par.outlets().size(); k++){
          rel_werror = std::abs((p_d_nonlinear[k] - p_d0_nonlinear[k])/p_d_nonlinear[k]);
        }
        MADprint("Navier Stokes: Non linear solver. Iteration = " + to_string(nn_iteration) + ". Windkessel residual = " + to_string(rel_werror) + "\n");

        // Update non-linear iteration
        VecZeroEntries(u_iter);
        VecAXPY(u_iter, 1.0, u);
        p_d0_nonlinear = p_d_nonlinear;

        if ((par.nonLinearTolerance() > rel_error) &&
            (par.nonLinearTolerance()/10.0 > rel_werror)){
          break;
        }
      }
      p_d0 = p_d0_nonlinear;
      VecDestroy(&u_iter);
      VecDestroy(&b_iter);

      // // Export windkessels pressure to txt
      // wkssl_vec.push_back(wkssl);
      // if (m_world_rank == 0){
      //   std::ofstream fout("wkssl_pressure.txt");
      //   fout << std::setprecision(10);
      //   for(auto const& x : wkssl_vec){
	    //     for (auto const& y : x){
      //       fout << y << " ";
      //     }
      //     fout << std::endl;
      //   }
      // }
    }

  private:
    vector<double> p_d0;
    vector<double> p_d0_nonlinear;
    vector<double> p_d_nonlinear;
    vector<double> wkssl;
    vector<vector<double>> wkssl_vec;
    Boundary bd;
    Geometry geo;
    Parameters par;
    int world_rank;
    Calculus m_calculus;
    IO m_io;
    NS ns;
    InnerProduct m_ip;
    int m_world_rank;
};

