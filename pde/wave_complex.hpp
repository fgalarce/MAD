/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2023,
    
     Felipe Galarce at INRIA/WIAS/PUCV

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

#ifndef MAD_WAVE_COMPLEX
#define MAD_WAVE_COMPLEX

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
#include <calculus.hpp>
#include <fem.hpp>

using namespace std;

class WaveComplex{

  public:
    WaveComplex(){}
    ~WaveComplex(){}

    void initialize(Parameters parameters, const Geometry & geometry, const Boundary & boundary);
    void finalize();
    void setInitialCondition();
    void setSolver();
    void setLHS(Mat A, Mat P = NULL);
    void solve(Vec b, Vec u);

    Mat assembleLHS();
    Vec assembleRHS();
    Mat A;
    KSP ksp;
    Calculus calculus;

    int nbDofs;
  private:
    Geometry geo;
    Boundary bd;
    IO io;
    MasterElement fe;
    MasterElementBD feBD;
    Parameters par;
    
    PetscErrorCode code;
    int m_verbose = 4;

    int nbDofsPerNode; 
    int nbVertices;
    int nbNodesPerElement;
    vector<int> nbDofVar;

    int m_world_rank, m, n;
    
    FEM fem;
};  

#endif
