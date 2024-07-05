/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2024,
    
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

#ifndef MAD_elasticity
#define MAD_elasticity 

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

using namespace std;

class Elasticity{

  public:
    Elasticity(){}
    ~Elasticity(){}

    void initialize(Parameters parameters, const Geometry & geometry);
    void finalize();
    Mat assembleLHS();
    Vec assembleRHS();
    void setSolver();
    void setLHS(Mat C);
    void solve(Vec b, Vec u);

    Mat A, M;
    KSP ksp;

  private:
    Geometry geo;
    MasterElement fe;
    Parameters par;

    PetscErrorCode code;

    int nbDofsPerNode; 
    int nbVertices;
    int nbDofs;

    int m_world_rank, m_verbose;


};  

#endif
