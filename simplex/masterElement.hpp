/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2021,
    
     Felipe Galarce at INRIA and WIAS

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

#ifndef MAD_Element
#define MAD_Element

#include <iostream>
#include <math.h>
#include <vector>
#include <slepc.h>
#include <STLvectorUtils.hpp>
#include <petscWrapper.hpp>
#include <tools.hpp>
#include <parameters.hpp>

using namespace std;

class MasterElement{

  public:
    MasterElement(){}
    ~MasterElement(){}

    void initialize(Parameters parameters, int dimension);
    void setCoordinates(const vector<vector<double>> & coordinates);
    void computeSize();
    void computeVolume();
    bool isInside(const vector<double> & point);
    double mass(int i, int j);
    double dphi_dx(int phi_id, int coord);
    double stiff(int i, int j);
    double stiffness(int i, int j);
    double stiffness_symmetric(int i, int j);
    double advection(int i, int j, vector<double> a);
    double mixed(int i, int j, int comp); // deprecated
    double dphi_i_phi_j(int i, int j, int comp);
    vector<double> dphi_i_phi_j(int i, int j);
    vector<double> dphi_i_phi_j_phi_i(int i, int j);
    double div_div(int i, int j, int comp);
    double volume();
    double phi(const vector<double> & point, int id_shape_function);
    double phi(const int & qp_point_id, int id_shape_function);
    vector<double> grad_phi(int id_shape_function);
    vector<double> barycentric_coor(vector<double> point);

    /* CFD */
    double a_grad_phi_i__dot__a_grad_phi_j(int i, int j, vector<double> a);

    /* Acces functions */
    inline const double & size() const {
      return m_size;}
    inline const vector<vector<double>> & jacobianInverse() const {
      return m_J_inv;}
    inline const double & detJacobian() const {
      return m_detJac;}
    inline const vector<double> & weights() const {
      return m_weights;}
    inline const double & weightsSum() const {
      return m_weights_sum;}
    inline const double & volume() const {
      return m_volume;}
    inline const vector<vector<double>> & grad_phi() const {
      return m_grad_phi;
    };

  private:
    Parameters par;
    /* Shape functions */
    vector<vector<double>> m_phi;
    vector<vector<double>> m_grad_phi; 
    vector<vector<double>> m_dphi_dx;
    vector<vector<double>> m_coordinates;
    double m_detJac;
    double m_size;
    double m_volume;
    vector<vector<double>> m_J; /* jacobian */
    vector<vector<double>> m_J_inv; /* jacobian inverse */
    vector<vector<double>> m_quadraturePoints;
    vector<double> m_weights;
    double m_weights_sum; /* Convenient to have when computing stiffness */
    vector<double> barycentric_coor2D(vector<double> point);
    vector<double> barycentric_coor3D(vector<double> point);
//    void arrangeReferenceDomain();

    int m_dimension;
    int m_nbNodesPerElement;
    int m_nbDofsPerNode;
    int m_nbDofsPerElement;

    PetscErrorCode code;
    int m_world_rank;

  public:
    vector<vector<double>> m_elementMass;
    double phi_i_phi_j(int i, int j);
    double dphi_i_COMP_dphi_j_COMP(int i, int j, int comp_i, int comp_j);
    double phi_i_dphi_j_COMP(int i, int j, int comp_j);
    double phi_dphi_i_COMP_phi_j(int i, int j, int x, vector<double> u);
    double phi_phi_dphi_i_COMP_dphi_j_COMP(int i, int j, int x, int y, vector<double> u, vector<double> v);
    double phi_phi_i_phi_j(int i, int j, vector<double> u);

};  

#endif
