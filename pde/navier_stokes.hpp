/*=============================================================================
  This file is part of the code MAD
  Copyright (C) 2017-2023,
    
     Felipe Galarce

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

#ifndef MAD_NS
#define MAD_NS

#include <iostream>
#include <math.h>
#include <vector>
#include <slepc.h>
#include <STLvectorUtils.hpp>
#include <petscWrapper.hpp>
#include <tools.hpp>
#include <parameters.hpp>
#include <masterElement.hpp>
#include <geometry.hpp>
#include <boundaries.hpp>
//#include <innerProduct.hpp>
#include <calculus.hpp>
#include <fem.hpp>

using namespace std;

class NavierStokes{

  public:
    NavierStokes(){}
    ~NavierStokes(){}

    void initialize(Parameters parameters, const Geometry & geometry, const Boundary & boundary);
    void finalize();
    void setInitialCondition();
    void update(Vec u, double time);
    void setSolver();
    void setLHS(Mat A, Mat P = NULL);
    void solve(Vec b, Vec u);
    void setInitialCondition(Vec up_0){
      u0 = zeros(nbDofs);
      VecAXPY(u0, 1.0, up_0);
    }
    inline const vector<double> & viscosity() const {
      return m_viscosity;
    }
    Vec viscosity_p1()  {
      return m_p1Viscosity;
    }

    double innerProduct(Vec u, Vec v);
    void assembleInnerProduct();

//    InnerProduct m_ip;
    Mat assembleLHS_static();
    Mat assembleLHS(Mat LHS_static, Vec u0_custom = NULL);
    Vec assembleRHS(Vec u0_custom = NULL);
    Mat A, M, C, K, B, C_rhs;
    KSP ksp;
    Calculus calculus;

    /* post proc */
    void computeFlows(Vec u);
    Vec shearStress(Vec u, int bdLabel);
    vector<double> getViscosity(Vec u0_custom = NULL);

    Vec u0;
    int nbDofs;
  private:
    Geometry geo;
    MasterElement fe;
    MasterElementBD feBD;
    Parameters par;
    Boundary bd;
    vector<double> m_viscosity;
    Vec m_p1Viscosity;
    INT m_interpolator;

    PetscErrorCode code;
    int m_verbose;

    int nbDofsPerNode; 
    int nbVertices;
    int nbDofsVel;
    int nbDofsPress;
    int nbNodesPerElement;

    bool m_non_newtonian = false;

    int m_world_rank, m, n;
    vector<Mat> m_ip;

    FEM femStat;
    FEM femMass;
    FEM femStiff;
    FEM femConv;
    FEM femNonSym;
    FEM femRHS;

};  

#endif
