/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2021,
    
     Felipe Galarce at INRIA

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

#ifndef MAD_laplacian
#define MAD_laplacian

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
#include <fem.hpp>
#include <boundaries.hpp>

using namespace std;

class Laplacian{

  public:
    Laplacian(){}
    ~Laplacian(){}

    void initialize(Parameters parameters, const Geometry & geometry, const Boundary & boundary);
    void finalize();
    Mat assembleLHS();
    void setSolver();
    void setLHS(Mat C);
    void solve(Vec b, Vec u);
    inline const Mat mass() const {
      return M;}

    Mat A, M;
    KSP ksp;

  private:
    Geometry geo;
    MasterElement fe;
    Parameters par;
    Boundary bd;

    PetscErrorCode code;

    int nbNodesPerElement;
    int nbDofsPerNode; 
    int nbVertices;
    int nbDofs;

    int m_world_rank, m_verbose;
    int m,n;

    FEM femStat;
    FEM femMass;
    void grad_u_grad_v(vector<int> simplex, double coef, int testLabel = 0, int trialLabel = 0);
};  

#endif
