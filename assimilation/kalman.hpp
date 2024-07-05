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

#ifndef ULTRA_4D_FLOW_KALMAN
#define ULTRA_4D_FLOW_KALMAN

#include <iostream>
#include <math.h>
#include <vector>

#include <slepc.h>
#include <STLvectorUtils.hpp>
#include <petscWrapper.hpp>
#include <tools.hpp>
#include <geometry.hpp>
#include <innerProduct.hpp>
#include <measures.hpp>
#include <model.hpp>

#include <parameters.hpp>

using namespace std;

class Kalman {

  public:
    Kalman(){}
    ~Kalman(){}

    void initialize(Parameters parameters, const Geometry & geometry, const InnerProduct & innerProduct, const Measures usMeasures);
    void initialize(Parameters parameters, const Geometry & geometry, Model model, const InnerProduct & innerProduct, const Measures usMeasures);
    void initialize(Parameters parameters, const Geometry & geometry);

//    inline const vector<double> & delta() const {
//      return m_delta[m_indexCurrentBasis];}

    Vec mean(const vector<Vec> & sample);
    Mat covariance(const vector<Vec> & sample);
    Mat propagateCovariance(Mat C, Mat A);
    Mat assembleFilter(Mat C10, Mat G, Mat Gt, Mat T);
    vector<Vec> samples(int iteration);
    vector<Vec> m_phi;
  private:
    Parameters par;

    /* Ambient space H */
    Geometry geo;
    int nbDofs;

    /* ROM V_n */
    Model Vn;
    int m_nbMeasures;
    int m_nbModes;


    InnerProduct m_ip;

};  

#endif


