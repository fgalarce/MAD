/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2023,
    
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

#ifndef CALCULUS 
#define CALCULUS 

#include <iostream>
#include <math.h>
#include <vector>
#include <slepc.h>
#include <STLvectorUtils.hpp>
#include <petscWrapper.hpp>
#include <tools.hpp>
#include <parameters.hpp>
#include <masterElement.hpp>
#include <io.hpp>
#include <geometry.hpp>
#include <innerProduct.hpp>
#include <interpolate.hpp>
#include <boundaries.hpp>

using namespace std;

class Calculus{

  public:
    Calculus(){}
    ~Calculus(){}
    void finalize();

    void initialize(Parameters parameters, const Geometry & geometry, const Boundary & boundary);

    /* compute derivatives */
    Vec divergence(Vec u);
    Vec gradient(Vec u, int partId = 0);
    Vec gradientOnBoundary(Vec u, int bdLabel);
    Vec wallShearStress(vector<Vec> u, int bdLabel);

    /* compute integrals */
    double boundaryIntegralScalar(Vec u, int bdLabel);

    /* tools */
    vector<Vec> split(Vec u);
    Vec magnitude(Vec u);
    Vec join(vector<Vec> u, IO & inputOutput); /* split^-1*/
    vector<Vec> decompose_vector_field(Vec u, string mode = "contiguous");
    Vec recompose_vector_field(vector<Vec> u, string mode = "contiguous"); /* decompose_vector_field^{-1} */

  private:

    Vec gradientScalar(Vec u, int partId);
    Vec gradientOnBoundaryForVectorField(Vec u, int bdLabel);
    Vec gradientOnBoundaryForScalarField(Vec u, int bdLabel);

    Parameters par;
    IO io; 
    Geometry geo;
    Boundary bd;
    INT m_interpolator;
    MasterElement fe;
    MasterElementBD feBD;

    int m_dimension;
    int m_nbDofs;
    int m_nbVertices;

    PetscErrorCode code;
    int m_verbose = 4;
    int m_world_rank;

};  

#endif
