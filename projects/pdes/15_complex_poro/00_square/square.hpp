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

vector<double> manufacturedURe(vector<double> x, double t, Parameters par){
  vector<double> bc(2);
  bc[0] = (x[0] - 1.0)*(x[0] - 1.0) * x[1]*x[1];
  bc[1] = x[0]*x[1]*(x[0] - 1.0);
  return bc;
}

vector<double> manufacturedUIm(vector<double> x, double t, Parameters par){
  vector<double> bc(2);
  bc[0] = (x[0] - 1.0) * (x[0] + 2.0) * (x[0] + 2.0) * x[1] * (x[1] + 1);
  bc[1] = 2.0 * x[0] * x[0] * x[1] * (x[0] - 1.0);  
  return bc;
}

vector<double> manufacturedPRe(vector<double> x, double t, Parameters par){
  vector<double> bc(1);
  bc[0] = sin(PI*x[0]/2.0) * cos(PI*x[1]/2.0);
  return bc;
}

vector<double> manufacturedPIm(vector<double> x, double t, Parameters par){
  vector<double> bc(1);
  bc[0] = (1.0 - cos(PI*x[0])) * (1.0 + cos(PI*x[1]));
  return bc;
}

vector<double> manufacturedPhiRe(vector<double> x, double t, Parameters par){
  vector<double> bc(1);
  int partId = 0;
  double lambda = par.youngModulus()[partId] * par.poissonRatio()[partId] / ((1.0 + par.poissonRatio()[partId])*(1.0 - 2.0*par.poissonRatio()[partId]));
  double trEpsUre = (x[0] - 1.0)*x[1]*x[1] + x[0] * (x[0] - 1.0);
  bc[0] = sin(PI*x[0]/2.0) * cos(PI*x[1]/2.0) - lambda * trEpsUre;
  return bc;
}

vector<double> manufacturedPhiIm(vector<double> x, double t, Parameters par){
  vector<double> bc(1);
  int partId = 0;
  double lambda = par.youngModulus()[partId] * par.poissonRatio()[partId] / ((1.0 + par.poissonRatio()[partId])*(1.0 - 2.0*par.poissonRatio()[partId]));
  double trEpsUim = x[0]*x[0]*(x[0] - 1.0) + x[1]*(x[0] - 1.0)*(x[0] + 2.0)*(x[1]+1.0);
  bc[0] = (1.0 - cos(PI*x[0])) * (1.0 + cos(PI*x[1])) - lambda * trEpsUim;
  return bc;
}


class Square{

  public:

    Square(){}
    ~Square(){};

    void initialize(Boundary & boundary, Geometry & geometry, Calculus calculus, Parameters parameters){
      par = parameters;
      bd = boundary;
      geo = geometry;
      cal = calculus;
    }

    /* update time, declare and impose boundary conditions */
    void applyBC(Mat A, Vec b){
      /* right and bottom for displacement  */
      bd.Dirichlet(par.walls()[0], manufacturedURe, 0);
      bd.Dirichlet(par.walls()[1], manufacturedURe, 0);
      bd.Dirichlet(par.walls()[0], manufacturedUIm, 1);
      bd.Dirichlet(par.walls()[1], manufacturedUIm, 1);
      /* left and top for pressure */
      bd.Dirichlet(par.walls()[2], manufacturedPRe, 2);
      bd.Dirichlet(par.walls()[3], manufacturedPRe, 2);
      bd.Dirichlet(par.walls()[2], manufacturedPIm, 3);
      bd.Dirichlet(par.walls()[3], manufacturedPIm, 3);
      bd.block(A);
      bd.block(b);
    }

    Vec manufacture(PoroComplex & poro){
//      Vec b = poro.computePseudoPressure(manufacturedURe, manufacturedUIm, manufacturedPRe, manufacturedPIm);
      Vec b = poro.assembleRHS(manufacturedPhiRe, 4);
      VecAXPY(b, 1.0, poro.assembleRHS(manufacturedPhiIm, 5));
      VecAXPY(b, 1.0, poro.assembleRHS(manufacturedURe, 0));
      VecAXPY(b, 1.0, poro.assembleRHS(manufacturedUIm, 1));
      VecAXPY(b, 1.0, poro.assembleRHS(manufacturedPRe, 2));
      VecAXPY(b, 1.0, poro.assembleRHS(manufacturedPIm, 3));
      
      code = VecAssemblyBegin(b); CHKERR(code);
      code = VecAssemblyEnd(b); CHKERR(code);
      return b;
    }

    double compare_anal(Vec u, IO & io){
      ip.initialize(par, geo, bd);
      int nbDofs = 0;
      nbDofVar.resize(par.nbVariables());
      for (int i = 0; i < par.nbVariables(); i++){
        nbDofs += par.nbDofsPerNode()[i]*geo.nbVertices;
        nbDofVar[i] = par.nbDofsPerNode()[i]*geo.nbVertices;
      }
      Vec sol_gt = vec(nbDofs);
      Vec u_gt = vec(nbDofs);
      assembleGT(u_gt, io);
      Vec error = zeros(nbDofs);
      VecAXPY(error, 1.0, u);
      VecAXPY(error, -1.0, u_gt);
//      return sqrt(ip(error, error));
      return norm(error);
    }

    double compare_anal(Vec u, PoroComplex & poro, IO & io){

      IO io_gt;
      io_gt.initialize(par.geometry(), "./results/");

      Geometry geo_gt;
      geo_gt.initialize(io_gt);

      INT interpolator;
      interpolator.initialize(par, geo_gt, io_gt);

      Boundary boundary;
      boundary.initialize(par, geo_gt);

      Calculus calculus;
      calculus.initialize(par, geo_gt, boundary);

      vector<Vec> sol_int(par.nbVariables());
      for (int i = 0; i < par.nbVariables(); i++){
        sol_int[i] = interpolator.interpolate_field(cal.split(u)[i], io, i); 
      }

      int nbDofsGT = 0;
      nbDofVar.resize(par.nbVariables());
      for (int i = 0; i < par.nbVariables(); i++){
        double normSolINT = norm(sol_int[i]);
        MADprint("normSOl_int:", normSolINT );
        nbDofsGT += par.nbDofsPerNode()[i]*geo_gt.nbVertices;
        nbDofVar[i] = par.nbDofsPerNode()[i]*geo_gt.nbVertices;
      }
      Vec sol_gt = vec(nbDofsGT);
      assembleGT(sol_gt, io_gt);

      io_gt.writeState(sol_int[0], "Iure");
      io_gt.writeState(sol_int[1], "Iuim");
      io_gt.writeState(sol_int[2], "Ipre");
      io_gt.writeState(sol_int[3], "Ipim");
      io_gt.writeState(sol_int[4], "Iphir");
      io_gt.writeState(sol_int[5], "Iphim");

      io_gt.writeState(calculus.split(sol_gt)[0], "GTu_re");
      io_gt.writeState(calculus.split(sol_gt)[1], "GTu_im");
      io_gt.writeState(calculus.split(sol_gt)[2], "GTp_re");
      io_gt.writeState(calculus.split(sol_gt)[3], "GTp_im");
      io_gt.writeState(calculus.split(sol_gt)[4], "GTphi_re");
      io_gt.writeState(calculus.split(sol_gt)[5], "GTphi_im");

      for (int i = 0; i < par.nbVariables(); i++){
        double normSolGT = norm(calculus.split(sol_gt)[i]);
        MADprint("normSOl_gt:", normSolGT );
      }

      /* Joint error */
      Vec sol_int_joint = cal.join(sol_int, io_gt);
      Vec error = zeros(nbDofsGT);
      VecAXPY(error, 1.0, sol_int_joint);
      VecAXPY(error, -1.0, sol_gt);
      ip.initialize(par, geo_gt, boundary);
      double normH1error = sqrt(ip(error, error));

      poro.assembleInnerProduct(geo_gt);
      vector<double> errorVar(nbDofVar.size());
      vector<Vec> errorSplit(nbDofVar.size());
      for (int i = 0; i < nbDofVar.size(); i++){
        vec(errorSplit[i], nbDofVar[i]);
        VecAXPY(errorSplit[i], 1.0, sol_int[i]);
        VecAXPY(errorSplit[i], -1.0, calculus.split(sol_gt)[i]);
        errorVar[i] = sqrt(poro.innerProduct(errorSplit[i], errorSplit[i], geo_gt)); 
        MADprint("Error in var " + par.variableName()[i] + ": ", errorVar[i]);
      }

      return normH1error;
    }

    void assembleGT(Vec u_gt, IO & io_gt){
      for (int i = 0; i < io_gt.nbVertices(); i++){
        vector<double> uRe_manu = manufacturedURe(io_gt.coordinates()[i], 0.0, par);
        vector<double> uIm_manu = manufacturedUIm(io_gt.coordinates()[i], 0.0, par);
        vector<double> pRe_manu = manufacturedPRe(io_gt.coordinates()[i], 0.0, par);
        vector<double> pIm_manu = manufacturedPIm(io_gt.coordinates()[i], 0.0, par);
        vector<double> phiRe_manu = manufacturedPhiRe(io_gt.coordinates()[i], 0.0, par);
        vector<double> phiIm_manu = manufacturedPhiIm(io_gt.coordinates()[i], 0.0, par);
        for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
          code = VecSetValue(u_gt, par.nbDofsPerNode()[0]*i+comp, uRe_manu[comp], INSERT_VALUES); CHKERR(code);
          code = VecSetValue(u_gt, nbDofVar[0] + par.nbDofsPerNode()[1]*i+comp, uIm_manu[comp], INSERT_VALUES); CHKERR(code);
        }
        code = VecSetValue(u_gt, nbDofVar[0] + nbDofVar[1] + i, pRe_manu[0], INSERT_VALUES); CHKERR(code);
        code = VecSetValue(u_gt, nbDofVar[0] + nbDofVar[1] + nbDofVar[2] + i, pIm_manu[0], INSERT_VALUES); CHKERR(code);
        code = VecSetValue(u_gt, nbDofVar[0] + nbDofVar[1] + nbDofVar[2] + nbDofVar[3] + i, phiRe_manu[0], INSERT_VALUES); CHKERR(code);
        code = VecSetValue(u_gt, nbDofVar[0] + nbDofVar[1] + nbDofVar[2] + nbDofVar[3] + nbDofVar[4] + i, phiIm_manu[0], INSERT_VALUES); CHKERR(code);
      }
      VecAssemblyBegin(u_gt);
      VecAssemblyEnd(u_gt);
    }

  private:
    Boundary bd;
    Parameters par;
    Geometry geo;
    Calculus cal;
    InnerProduct ip;
    PetscErrorCode code;
    vector<int> nbDofVar;
};
