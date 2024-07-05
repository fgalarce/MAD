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

#ifndef MAD_ROM
#define MAD_ROM

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <random>
#include <chrono>
#include <slepc.h>
#include <STLvectorUtils.hpp>

#include <petscWrapper.hpp>
#include <tools.hpp>
#include <parameters.hpp>
#include <geometry.hpp>
#include <io.hpp>
#include <innerProduct.hpp>
#include <boundaries.hpp>
#include <linearAlgebra.hpp>
#include <calculus.hpp>

using namespace std;

class ROM {

  public:
    ROM(){}
    ~ROM(){}

    void initialize(Parameters parameters, Geometry & geometry, IO & inputOutput, Calculus & calc);
    void finalize();

    inline const Vec basis(int i) const {
      return m_basis[i];}

  private:

    Parameters par;
    PetscErrorCode code;

    int m_nbDofs;
    int m_world_rank;
    vector<Vec> m_basis;

    Geometry geo;
    IO io;
    Calculus calculus;

};  

#endif
