/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2020,
    
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

#ifndef MAD_Interpolate
#define MAD_Interpolate

#include <iostream>
#include <math.h>
#include <vector>
#include <slepc.h>
#include <STLvectorUtils.hpp>
#include <petscWrapper.hpp>
#include <tools.hpp>
#include <parameters.hpp>
#include <masterElementBD.hpp>
#include <masterElement.hpp>
#include <io.hpp>
#include <geometry.hpp>

using namespace std;

class INT{

  public:
    INT(){}
    ~INT(){}

    void initialize(Parameters parameters, const Geometry & geometry);
    void initialize(Parameters parameters, const Geometry & geometry, const IO & inputOutput);
    vector<Vec> decompose_vector_field(Vec u);

    /* Field interpolation */
    vector<Vec> interpolate_field(vector<Vec> u_i, IO & io_i, int varLabel = 0);
           Vec  interpolate_field(       Vec  u_i, IO & io_i, int varLabel = 0);
           Vec  interpolate_field(string filePath, IO & io_i, int varLabel = 0);

    /* P0 - P1 interpolation routines */
    Vec interpolateP0_P1(vector<double> u);
    Vec interpolateP0_P1_boundary(vector<double> u, string mode = "splitted");

  private:

    Geometry geo;
    IO io;
    Parameters par;

    int m_dimension;
    int m_nbDofs;
    int m_nbVertices;

    vector<int> m_indexes; 
    vector<int> m_indexesBoundary; 

    PetscErrorCode code;
    int m_verbose = 10;
    int m_world_rank;

    MasterElement fe;
    MasterElementBD feBD;
};  


#endif
