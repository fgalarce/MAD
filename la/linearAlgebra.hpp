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


#ifndef ULTRA_4D_FLOW_LINEAR_ALGEBRA
#define ULTRA_4D_FLOW_LINEAR_ALGEBRA

#include <slepc.h>
#include <petscWrapper.hpp>
#include <parameters.hpp>
#include <geometry.hpp>
#include <innerProduct.hpp>

using namespace std;

class LinearAlgebra{

  public:
    LinearAlgebra(){}
    ~LinearAlgebra(){}
    void initialize(Parameters parameters, InnerProduct iproduct);
    /* If not InnerProduct is passed, l2 space is assumed */
    void initialize(Parameters parameters);
    bool checkOrthonormality(vector<Vec> u);
    vector<double> orthonormalize(vector<Vec> & u); /* it returns an array with the norms */
    vector<vector<double>> eigenValueProblem(Mat & A, Mat & B, string mode);
    vector<vector<double>> eigenValueProblem(Mat & A, string mode="HEP");
    Vec projectOn(Vec u, vector<Vec> v);
    Mat & pseudoinverse(Mat A);

    InnerProduct ip;
  private:
    Parameters par;

    int delta(unsigned int i, unsigned int j){
      if (i == j){
        return 1;
      } else {
        return 0;
      }
    }

    PetscErrorCode code;
    int m_world_rank;
};  

#endif
