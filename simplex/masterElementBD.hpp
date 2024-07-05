/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2021,
    
     Felipe Galarce at INRIA / WIAS

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

#ifndef MAD_ElementBD
#define MAD_ElementBD

#include <iostream>
#include <math.h>
#include <vector>
#include <slepc.h>
#include <STLvectorUtils.hpp>
#include <petscWrapper.hpp>
#include <tools.hpp>
#include <parameters.hpp>

using namespace std;

class MasterElementBD{

  public:
    MasterElementBD(){}
    ~MasterElementBD(){}

    void initialize(Parameters parameters, int dimension);

    double phi(const vector<double> & point, int id_shape_function);
    double phi(const int & qp_point_id, int id_shape_function);
    vector<double> grad_phi(int id_shape_function){
      return m_grad_phi[id_shape_function];
    }

    double mass(int i, int j);
    double stiffness(int i, int j);
    double stiffness_symmetric(int i, int j);
    double mixed(int i, int j, int comp);
    double advection(int i, int j, vector<double> a);
    
    void setCoordinates(const vector<vector<double>> & coordinates);
    void computeSize();
    void computeVolume();

    /* Acces functions */
    inline const vector<vector<double>> & quadraturePoints() const {
      return m_quadraturePoints;}
    inline const vector<double> & weights() const {
      return m_weights;}
    inline const double & size() const {
      return m_size;}
    inline const double & volume() const {
      return m_volume;}
    inline const vector<vector<double>> & jacobianInverse() const {
      return m_J_inv;}
    void computeNormal();

    inline const vector<double> & normal() const {
      return m_normal;}

    double computeFlow(vector<double> vectorField);

    inline const double & sqrtDetJTJ() const { 
      return m_sqrtDetJTJ;
    }
  private:

    vector<vector<double>> m_quadraturePoints;
    vector<double> m_weights;
    double m_weights_sum;
//    void arrangeReferenceDomain();

    vector<vector<double>> m_grad_phi; 
    vector<vector<double>> m_J;
    vector<vector<double>> m_J_inv;
    vector<vector<double>> m_J_adj;
    vector<vector<double>> m_JTJ;
    vector<double> m_normal;
    double m_volume;

    Parameters par;

    int m_nbDofsPerNode;
    int m_nbNodesPerElement;
    int m_nbDofsPerElement;
    int m_dimension;
    double m_sqrtDetJTJ;

    vector<vector<double>> m_phi;
    vector<vector<double>> m_coordinates;
    double m_size;

    int m_world_rank;
};  

#endif
