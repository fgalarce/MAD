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

#ifndef MAD_FEM_
#define MAD_FEM_

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

using namespace std;

class FEM{

  public:
    FEM(){}
    ~FEM(){}

    void initialize(Parameters parameters, const Geometry & geometry, const Boundary & boundary, Mat & A);
    void finalize();

    void setSimplex(vector<int> element);
    void copySimplex(FEM & finiteElement);
    void du_x_dv_y                (double coef, int x, int y, int testLabel = 0, int trialLabel = 0);
    void u_dv_y                   (double coef, int y, int testLabel = 0, int trialLabel = 0);
    void du_x_v                   (double coef, int x, int testLabel = 0, int trialLabel = 0);
    void a_du_x_v                 (double coef, int x, vector<double> a, int testLabel, int trialLabel);
    void a_du_x_b_dv_y            (double coef, int x, int y, vector<double> a, vector<double> b, int testLabel, int trialLabel);

    void epsilon_u_epsilon_v      (double coef, int testLabel = 0, int trialLabel = 0);
    void u_dot_v                  (double coef, int testLabel = 0, int trialLabel = 0);
    void grad_u_grad_v            (double coef, int testLabel = 0, int trialLabel = 0);
    void u_dot_v_scalar           (double coef, int testLabel = 0, int trialLabel = 0);
    void u_div_v                  (double coef, int testLabel, int trialLabel);
    void div_u_v                  (double coef, int testLabel, int trialLabel);
    void div_u_div_v              (double coef, int testLabel, int trialLabel);
    void u__dot__grad_q           (double coef, int testLabel, int trialLabel);
    void a_grad_u_dot_v           (double coef, vector<double> u_el, int testLabel, int trialLabel);
    void a_grad_u__dot__a_grad_v  (double coef, vector<double> u_el, int testLabel, int trialLabel);
    void a_grad_v__dot__grad_p    (double coef, vector<double> u_el, int testLabel, int trialLabel);
    void a_grad_u__dot__grad_q    (double coef, vector<double> u_el, int testLabel, int trialLabel);
    void u__dot__a_grad_v         (double coef, vector<double> u_el, int testLabel, int trialLabel);
    inline double feSize() const {
      return fe.size(); }
    MasterElement finiteElement() {
      return fe; }
    vector<int> simplex(){
      return m_simplex; 
    }
    Mat massMatrix(){
      return A;
    }

  private:

    Mat A;
    int nbDofs, nbVertices, nbNodesPerElement;
    vector<int> nbDofVar;

    Geometry geo;
    Boundary bd;

    MasterElement fe;
    MasterElementBD feBD;
    Parameters par;

    PetscErrorCode code;
    int m_verbose;

    int m_world_rank, m, n;
    vector<int> m_simplex;
    
};  

#endif
