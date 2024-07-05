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

#ifndef MAD_PORO_COMPLEX
#define MAD_PORO_COMPLEX

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

class PoroComplex{

  public:
    PoroComplex(){}
    ~PoroComplex(){}

    void initialize(Parameters parameters, const Geometry & geometry, const Boundary & boundary, const IO & inputOutput);
    void finalize();
    void setInitialCondition();
    void setSolver();
    void setLHS(Mat A, Mat P = NULL);
    void solve(Vec b, Vec u);

    double innerProduct(Vec u, Vec v, Geometry & geometry);
    void assembleInnerProduct(Geometry & geometry);
    inline const Vec analytic_solution() const {
      return m_uAnal;}

    void manufacture(Vec sol, int offset, vector<double> (*manufacture) (vector<double>,double,Parameters), int idVariable);
    Vec computePseudoPressure(vector<double> (*Ure) (vector<double>,double,Parameters),
                              vector<double> (*Uim) (vector<double>,double,Parameters), 
                              vector<double> (*Pre) (vector<double>,double,Parameters),
                              vector<double> (*Pim) (vector<double>,double,Parameters));
    Mat assembleLHS();
    Vec assembleRHS(vector<double> (*manufactured_solution) (vector<double>,double,Parameters), int idVariable);
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
    
    vector<double> m_viscosity;
    Vec m_p1Viscosity;
    INT m_interpolator;

    PetscErrorCode code;
    int m_verbose = 4;

    int nbDofsPerNode; 
    int nbVertices;
    int nbNodesPerElement;
    vector<int> nbDofVar;

    bool m_non_newtonian = false;

    int m_world_rank, m, n;
    Vec m_uAnal;
    Vec phi_re, phi_im;

    /* Ambient space ip */
    vector<Mat> m_ip;
    
    FEM fem;
};  

#endif
