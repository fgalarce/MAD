/*=============================================================================
  This file is part of the code MAD
  Copyright (C) 2017-2023,
    
     Felipe Galarce at INRIA / WIAS / PUCV

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

#define EPSILON 0.000000000001

#include<masterElement.hpp>
#include<numeric>

void MasterElement::initialize(Parameters parameters, int dimension){

  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0 and par.debug() == true) cout << "Element: Initializing." << endl;

  par = parameters;
  m_dimension = dimension;
  m_nbNodesPerElement = m_dimension + 1;
  m_nbDofsPerNode = par.nbDofsPerNode()[0];
  m_nbDofsPerElement = m_nbNodesPerElement*m_nbDofsPerNode;

  /*  Quadrature data for tetrahedron
      Refs P Keast, Moderate degree tetrahedral quadrature formulas, CMAME 55: 339-348 (1986)
      O. C. Zienkiewicz, The Finite Element Method,  Sixth Edition  
  */
  assert(((par.nbQuadraturePoints() == 3) or 
          (par.nbQuadraturePoints() == 4) or 
          (par.nbQuadraturePoints() == 5) or 
          (par.nbQuadraturePoints() == 1)) and 
          "Gauss quadrature rule not valid.");

  m_quadraturePoints.resize(par.nbQuadraturePoints());
  m_weights.resize(par.nbQuadraturePoints());
  for (int i = 0; i < par.nbQuadraturePoints(); i++){ 
    m_quadraturePoints[i].resize(m_dimension);  
  }

  if (m_dimension == 1){
    if (par.nbQuadraturePoints() == 1){
      m_quadraturePoints[0][0] = 0.5;
      m_quadraturePoints[0][1] = 0.0;
      m_quadraturePoints[0][2] = 0.0;
      m_weights[0] = 1.0;
    } else if (par.nbQuadraturePoints() == 3){
      m_quadraturePoints[0][0] = 0.5* 0 + 0.5;
      m_quadraturePoints[0][1] = 0.0;
      m_quadraturePoints[0][2] = 0.0;
      m_weights[0] = 4.0/3.0 * 0.5;

      m_quadraturePoints[1][0] = 0.5* 1 + 0.5;
      m_quadraturePoints[1][1] = 0.0;
      m_quadraturePoints[1][2] = 0.0;
      m_weights[1] = 1.0/3.0 * 0.5;

      m_quadraturePoints[2][0] = 0.5* -1 + 0.5;
      m_quadraturePoints[2][1] = 0.0;
      m_quadraturePoints[2][2] = 0.0;
      m_weights[2] = 1.0/3.0 * 0.5;
    } else {
      exit(1);
    }
  }
  if (m_dimension == 2){
    if (par.nbQuadraturePoints() == 1){
      // Scheme from Zienkiewicz and Taylor, 1 point, degree of precision 1
      m_quadraturePoints[0][0] = 0.333333333333333333;
      m_quadraturePoints[0][1] = 0.333333333333333333;
      m_quadraturePoints[0][2] = 0.0;
      m_weights[0] = 0.5;
    } else if(par.nbQuadraturePoints() == 3){
      // Scheme from Strang and Fix, 3 points, degree of precision 2
      m_quadraturePoints[0][0] = 1.0/6.0;
      m_quadraturePoints[0][1] = 1.0/6.0;
      m_quadraturePoints[0][2] = 0.0;
      m_weights[0] = 1.0/6.0;

      m_quadraturePoints[1][0] = 1.0/6.0;
      m_quadraturePoints[1][1] = 2.0/3.0;
      m_quadraturePoints[1][2] = 0.0;
      m_weights[1] = 1.0/6.0;

      m_quadraturePoints[2][0] = 2.0/3.0;
      m_quadraturePoints[2][1] = 1.0/6.0;
      m_quadraturePoints[2][2] = 0.0;
      m_weights[2] = 1.0/6.0;
    } else if(par.nbQuadraturePoints() == 4){
      m_quadraturePoints[0][0] = 1.0/3.0;
      m_quadraturePoints[0][1] = 1.0/3.0;
      m_quadraturePoints[0][2] = 0.0;
      m_weights[0] = -27.0/96.0;

      m_quadraturePoints[1][0] = 0.6;
      m_quadraturePoints[1][1] = 0.2;
      m_quadraturePoints[1][2] = 0;
      m_weights[1] = 25.0/96.0;

      m_quadraturePoints[2][0] = 0.2;
      m_quadraturePoints[2][1] = 0.6;
      m_quadraturePoints[2][2] = 0.0;
      m_weights[2] = 25.0/96.0;

      m_quadraturePoints[3][0] = 0.2;
      m_quadraturePoints[3][1] = 0.2;
      m_quadraturePoints[3][2] = 0.0;
      m_weights[3] = 25.0/96.0;
    } else if(par.nbQuadraturePoints() == 6){
      // Scheme from Strang and Fix, 6 points, degree of precision 3
      m_quadraturePoints[0][0] = 0.659027622374092;
      m_quadraturePoints[0][1] = 0.231933368553031;
      m_quadraturePoints[0][2] = 0.0;
      m_weights[0] = 1.0/12.0;

      m_quadraturePoints[1][0] = 0.659027622374092;
      m_quadraturePoints[1][1] = 0.109039009072877;
      m_quadraturePoints[1][2] = 0;
      m_weights[1] = 1.0/12.0;

      m_quadraturePoints[2][0] = 0.231933368553031;
      m_quadraturePoints[2][1] = 0.659027622374092;
      m_quadraturePoints[2][2] = 0.0;
      m_weights[2] = 1.0/12.0;

      m_quadraturePoints[3][0] = 0.231933368553031;
      m_quadraturePoints[3][1] = 0.109039009072877;
      m_quadraturePoints[3][2] = 0.0;
      m_weights[3] = 1.0/12.0;

      m_quadraturePoints[4][0] = 0.109039009072877;
      m_quadraturePoints[4][1] = 0.659027622374092;
      m_quadraturePoints[4][2] = 0.0;
      m_weights[4] = 1.0/12.0;

      m_quadraturePoints[5][0] = 0.109039009072877;
      m_quadraturePoints[5][1] = 0.231933368553031;
      m_quadraturePoints[5][2] = 0.0;
      m_weights[5] = 1.0/12.0;
    } else {
      exit(1);
    }
  }
  if (m_dimension == 3){
    if (par.nbQuadraturePoints() == 1){
      // Scheme from Zienkiewicz and Taylor, 1 point, degree of precision 1
      m_quadraturePoints[0][0] = 0.25;
      m_quadraturePoints[0][1] = 0.25;
      m_quadraturePoints[0][2] = 0.25;
      m_weights[0] = 1.0/6.0;

    } else if (par.nbQuadraturePoints() == 4){
      // Scheme from Zienkiewicz and Taylor, 4 points, degree of precision 2
      double a, b;
      a = 0.585410196624969;
      b = 0.138196601125011;
      m_quadraturePoints[0][0] = a;
      m_quadraturePoints[0][1] = b;
      m_quadraturePoints[0][2] = b;
      m_weights[0] = 1.0/24;

      m_quadraturePoints[1][0] = b;
      m_quadraturePoints[1][1] = a;
      m_quadraturePoints[1][2] = b;
      m_weights[1] = 1.0/24;

      m_quadraturePoints[2][0] = b;
      m_quadraturePoints[2][1] = b;
      m_quadraturePoints[2][2] = a;
      m_weights[2] = 1.0/24;

      m_quadraturePoints[3][0] = b;
      m_quadraturePoints[3][1] = b;
      m_quadraturePoints[3][2] = b;
      m_weights[3] = 1.0/24;
    } else if (par.nbQuadraturePoints() == 5){
      // # Scheme from Zienkiewicz and Taylor, 5 points, degree of precision 3
      // Note: this scheme has a negative weight
      m_quadraturePoints[0][0] = 0.25;
      m_quadraturePoints[0][1] = 0.25;
      m_quadraturePoints[0][2] = 0.25;
      m_weights[0] = -0.8/6.0;

      m_quadraturePoints[1][0] = 0.5;
      m_quadraturePoints[1][1] = 0.1666666666666667;
      m_quadraturePoints[1][2] = 0.1666666666666667;
      m_weights[1] = 0.45/6.0;

      m_quadraturePoints[2][0] = 0.1666666666666667;
      m_quadraturePoints[2][1] = 0.1666666666666667;
      m_quadraturePoints[2][2] = 0.1666666666666667;
      m_weights[2] = 0.45/6.0;

      m_quadraturePoints[3][0] = 0.1666666666666667;
      m_quadraturePoints[3][1] = 0.1666666666666667;
      m_quadraturePoints[3][2] = 0.5;
      m_weights[3] = 0.45/6.0;

      m_quadraturePoints[4][0] = 0.1666666666666667;
      m_quadraturePoints[4][1] = 0.5;
      m_quadraturePoints[4][2] = 0.1666666666666667;
      m_weights[4] = 0.45/6.0;
    } else {
      exit(1);
    }
  }

  /* Pre-compute shape functions on quadrature points */
  m_phi.resize(m_nbNodesPerElement);
  for (int i = 0; i < m_nbNodesPerElement; i++){
    m_phi[i].resize(par.nbQuadraturePoints());
  }
  for (int i = 0; i < par.nbQuadraturePoints(); i++){
    if (m_dimension == 3){
      m_phi[0][i] = 1 - m_quadraturePoints[i][0] - m_quadraturePoints[i][1] - m_quadraturePoints[i][2]; /* 1 - \xi - \eta - \gamma */
      m_phi[1][i] = m_quadraturePoints[i][0]; /* \xi    */
      m_phi[2][i] = m_quadraturePoints[i][1]; /* \eta   */
      m_phi[3][i] = m_quadraturePoints[i][2]; /* \gamma */
    } else if (m_dimension == 2){
      m_phi[0][i] = 1 - m_quadraturePoints[i][0] - m_quadraturePoints[i][1]; /* 1 - \xi - \eta  */
      m_phi[1][i] = m_quadraturePoints[i][0]; /* \xi    */
      m_phi[2][i] = m_quadraturePoints[i][1]; /* \eta   */
    } else if (m_dimension == 1){
      m_phi[0][i] = 1 - m_quadraturePoints[i][0]; /* 1 - \xi   */
      m_phi[1][i] = m_quadraturePoints[i][0]; /* \xi    */
    } else {
      exit(1);
    }
  } 

  /* Initialize shape functions derivatives, constants for P1 elements */
  m_grad_phi.resize(m_nbNodesPerElement);
  for (int j = 0; j < m_nbNodesPerElement; j++){
    m_grad_phi[j].resize(m_dimension);
  }

  if (m_dimension == 1){
    m_grad_phi[0][0] = -1;
    m_grad_phi[1][0] = 1;
  } else if (m_dimension == 2){
    m_grad_phi[0][0] = -1;
    m_grad_phi[0][1] = -1;

    m_grad_phi[1][0] = 1;
    m_grad_phi[1][1] = 0;

    m_grad_phi[2][0] = 0;
    m_grad_phi[2][1] = 1;
  } else if (m_dimension == 3){
    m_grad_phi[0][0] = -1;
    m_grad_phi[0][1] = -1;
    m_grad_phi[0][2] = -1;

    m_grad_phi[1][0] = 1;
    m_grad_phi[1][1] = 0;
    m_grad_phi[1][2] = 0;

    m_grad_phi[2][0] = 0;
    m_grad_phi[2][1] = 1;
    m_grad_phi[2][2] = 0;

    m_grad_phi[3][0] = 0;
    m_grad_phi[3][1] = 0;
    m_grad_phi[3][2] = 1;
//    m_grad_phi = transpose(m_grad_phi);
  }

  m_J.resize(m_dimension); 
  for (int i = 0; i < m_dimension; i++){
    m_J[i].resize(m_dimension);
  }

  m_weights_sum = 0.0;
  for (int qp = 0; qp < par.nbQuadraturePoints(); qp++){
    m_weights_sum += m_weights[qp];
  }

  m_elementMass.resize(m_nbNodesPerElement);
  for (int i=0; i < m_nbNodesPerElement; i++){
    m_elementMass[i].resize(m_nbNodesPerElement);
    for (int j=0; j < m_nbNodesPerElement; j++){
      m_elementMass[i][j] = 0.0;
      for (int qp = 0; qp < par.nbQuadraturePoints(); qp++){
        m_elementMass[i][j] += m_phi[i][qp] * m_phi[j][qp] * m_weights[qp];
      }
    }
  }
}

double MasterElement::phi_i_phi_j(int i, int j){
  return m_elementMass[i][j] * fabs(m_detJac);
}

double MasterElement::phi_phi_i_phi_j(int i, int j, vector<double> u){
  double conv = 0;
  for (int qp = 0; qp < par.nbQuadraturePoints(); qp++){
    for (int k = 0; k < m_dimension+1; k++){
      conv += u[k] * m_phi[k][qp] * m_phi[j][qp] * m_phi[i][qp] * m_weights[qp];
    }
  }
  return conv * fabs(m_detJac);
}


double MasterElement::phi_dphi_i_COMP_phi_j(int i, int j, int x, vector<double> u){
  double conv = 0;
  for (int qp = 0; qp < par.nbQuadraturePoints(); qp++){
    for (int k = 0; k < m_dimension+1; k++){
      conv += u[k] * m_phi[k][qp] * m_phi[j][qp] * m_weights[qp];
    }
  }
  return conv * m_dphi_dx[i][x] * fabs(m_detJac);
}


//double MasterElement::advection(int i, int j, vector<double> a){
//  double advection = 0;
//  for (int qp = 0; qp < par.nbQuadraturePoints(); qp++){
//    advection += phi(qp,i) * m_weights[qp];
//  }
//  return advection * dot(transpose(m_J_inv) * grad_phi(j), a) * m_detJac;
//}

double MasterElement::phi_phi_dphi_i_COMP_dphi_j_COMP(int i, int j, int x, int y, vector<double> u, vector<double> v){
  double conv = 0.0;
  for (int qp = 0; qp < par.nbQuadraturePoints(); qp++){
    for (int k = 0; k < m_dimension+1; k++){
      conv += u[k] * m_phi[k][qp] * v[k] * m_phi[k][qp] * m_weights[qp];
    }
  }
  return conv * m_dphi_dx[i][x] * m_dphi_dx[j][y] * fabs(m_detJac);
}

double MasterElement::dphi_i_COMP_dphi_j_COMP(int i, int j, int comp_i, int comp_j){
  return m_dphi_dx[i][comp_i] * m_dphi_dx[j][comp_j] * fabs(m_detJac) * m_weights_sum;
}

double MasterElement::phi_i_dphi_j_COMP(int i, int j, int comp_j){
  double mixed = 0.0;
  for (int qp = 0; qp < par.nbQuadraturePoints(); qp++){
    mixed += m_phi[i][qp] * m_weights[qp];
  }
  return m_dphi_dx[j][comp_j] * mixed * fabs(m_detJac);
}

void MasterElement::setCoordinates(const vector<vector<double>> & coordinates){

  m_coordinates = coordinates;
  if (m_dimension == 3){
    m_J[0][0] = m_coordinates[1][0] - m_coordinates[0][0];
    m_J[1][0] = m_coordinates[1][1] - m_coordinates[0][1];
    m_J[2][0] = m_coordinates[1][2] - m_coordinates[0][2];

    m_J[0][1] = m_coordinates[2][0] - m_coordinates[0][0];
    m_J[1][1] = m_coordinates[2][1] - m_coordinates[0][1];
    m_J[2][1] = m_coordinates[2][2] - m_coordinates[0][2];

    m_J[0][2] = m_coordinates[3][0] - m_coordinates[0][0];
    m_J[1][2] = m_coordinates[3][1] - m_coordinates[0][1];
    m_J[2][2] = m_coordinates[3][2] - m_coordinates[0][2];

  } else if (m_dimension == 2) {
    m_J[0][0] = m_coordinates[1][0] - m_coordinates[0][0]; 
    m_J[1][0] = m_coordinates[1][1] - m_coordinates[0][1];

    m_J[0][1] = m_coordinates[2][0] - m_coordinates[0][0];
    m_J[1][1] = m_coordinates[2][1] - m_coordinates[0][1];
  } else {
    m_J[0][0] = m_coordinates[1][0] - m_coordinates[0][0]; 
  }

  m_detJac = fabs(determinant(m_J));
  m_J_inv = invert(m_J);

  m_dphi_dx.resize(m_nbNodesPerElement);
  for (int i = 0; i < m_nbNodesPerElement; i++){
    m_dphi_dx[i].resize(m_dimension);
    m_dphi_dx[i] = transpose(m_J_inv) * m_grad_phi[i];
  }
}


double MasterElement::mass(int i, int j){
  double mass = 0.0;
  for (int qp = 0; qp < par.nbQuadraturePoints(); qp++){
    mass += m_phi[i][qp] * m_phi[j][qp] * m_weights[qp];
  }
  return mass * fabs(m_detJac);
}

double MasterElement::dphi_dx(int phi_id, int coord){
  return m_dphi_dx[phi_id][coord];
}

double MasterElement::stiff(int i, int j){
  return dot(m_dphi_dx[i], m_dphi_dx[j]) * fabs(m_detJac) * m_weights_sum;
}

double MasterElement::stiffness(int i, int j){
  return dot(transpose(m_J_inv) * m_grad_phi[i], transpose(m_J_inv) * m_grad_phi[j]) * fabs(m_detJac) * m_weights_sum;
}

double MasterElement::stiffness_symmetric(int i, int j){
  return 0.5 * (dot(transpose(m_J_inv) * m_grad_phi[i], transpose(m_J_inv) * m_grad_phi[j])
        + 2.0 * dot(transpose(transpose(m_J_inv) * m_grad_phi[i]), transpose(m_grad_phi[j]) * m_J_inv)
        + dot(transpose(m_grad_phi[i]) * m_J_inv, transpose(m_grad_phi[j]) * m_J_inv)) * fabs(m_detJac) * m_weights_sum;
}

double MasterElement::mixed(int i, int j, int comp){
  double mixed = 0.0;
  for (int qp = 0; qp < par.nbQuadraturePoints(); qp++){
    mixed += m_phi[j][qp] * m_weights[qp];
  }
  return m_dphi_dx[i][comp] * mixed * fabs(m_detJac);
}

double MasterElement::dphi_i_phi_j(int i, int j, int comp){
  double mixed = 0.0;
  for (int qp = 0; qp < par.nbQuadraturePoints(); qp++){
    mixed += m_phi[j][qp] * m_weights[qp];
  }
  return m_dphi_dx[i][comp] * mixed * fabs(m_detJac);
}

vector<double> MasterElement::dphi_i_phi_j(int i, int j){
  double mixed = 0.0;
  for (int qp = 0; qp < par.nbQuadraturePoints(); qp++){
    mixed += m_phi[j][qp] * m_weights[qp];
  }
  return m_dphi_dx[i] * mixed * fabs(m_detJac);
}

vector<double> MasterElement::dphi_i_phi_j_phi_i(int i, int j){
  double mixed = 0.0;
  for (int qp = 0; qp < par.nbQuadraturePoints(); qp++){
    mixed += m_phi[i][qp] * m_phi[j][qp] * m_weights[qp];
  }
  return m_dphi_dx[i] * mixed * fabs(m_detJac);
}

double MasterElement::div_div(int i, int j, int comp){
  double div_div = 0.0;
  for (int coord = 0; coord < m_nbDofsPerNode; coord++){
    div_div = div_div + (transpose(m_J_inv) * m_grad_phi[i])[coord];
  }
  return div_div*(transpose(m_J_inv) * m_grad_phi[j])[comp]*fabs(m_detJac);
}

double MasterElement::advection(int i, int j, vector<double> a){
  double advection = 0;
  for (int qp = 0; qp < par.nbQuadraturePoints(); qp++){
    advection += phi(qp,i) * m_weights[qp];
  }
  return advection * dot(transpose(m_J_inv) * grad_phi(j), a) * m_detJac;
}

vector<double> MasterElement::grad_phi(int id_shape_function){
  assert(id_shape_function < m_nbNodesPerElement);
  return m_grad_phi[id_shape_function];
}

double MasterElement::phi(const int & qp_point_id, int id_shape_function){
  return m_phi[id_shape_function][qp_point_id];
}

double MasterElement::volume(){
  if (m_dimension == 3){
    return 1.0/6.0 * abs(dot(m_coordinates[0] - m_coordinates[3], cross((m_coordinates[1] - m_coordinates[3]), (m_coordinates[2] - m_coordinates[3]))));
  } else{
    cout << "Implement volume for 2d" << endl;
    exit(1);
  }
}

void MasterElement::computeSize(){
  m_size = 0.0;
  for (int i = 0; i < m_nbNodesPerElement; i++){
    for (int j = 0; j < m_nbNodesPerElement; j++){
      double dist = norm(m_coordinates[i] - m_coordinates [j]);
      if (i != j && dist > m_size)
        m_size = dist;
    }
  }
}

void MasterElement::computeVolume(){
  if (m_dimension == 2){
    m_volume = 0.5 * norm(cross(m_coordinates[0] - m_coordinates[1], m_coordinates[0] - m_coordinates[2]));
  } else {
    m_volume = 1.0/6.0 * abs(dot(m_coordinates[0] - m_coordinates[3], cross((m_coordinates[1] - m_coordinates[3]), (m_coordinates[2] - m_coordinates[3]))));
  }
}

bool MasterElement::isInside(const vector<double> & point){

  vector<vector<double>> D0(4);
  D0[0] = m_coordinates[0];
  D0[0].push_back(1);
  D0[1] = m_coordinates[1];
  D0[1].push_back(1);
  D0[2] = m_coordinates[2];
  D0[2].push_back(1);
  D0[3] = m_coordinates[3];
  D0[3].push_back(1);

  vector<vector<double>> D1(4);
  D1[0] = point;
  D1[0].push_back(1);
  D1[1] = m_coordinates[1];
  D1[1].push_back(1);
  D1[2] = m_coordinates[2];
  D1[2].push_back(1);
  D1[3] = m_coordinates[3];
  D1[3].push_back(1);

  vector<vector<double>> D2(4);
  D2[0] = m_coordinates[0];
  D2[0].push_back(1);
  D2[1] = point;
  D2[1].push_back(1);
  D2[2] = m_coordinates[2];
  D2[2].push_back(1);
  D2[3] = m_coordinates[3];
  D2[3].push_back(1);

  vector<vector<double>> D3(4);
  D3[0] = m_coordinates[0];
  D3[0].push_back(1);
  D3[1] = m_coordinates[1];
  D3[1].push_back(1);
  D3[2] = point;
  D3[2].push_back(1);
  D3[3] = m_coordinates[3];
  D3[3].push_back(1);

  vector<vector<double>> D4(4);
  D4[0] = m_coordinates[0];
  D4[0].push_back(1);
  D4[1] = m_coordinates[1];
  D4[1].push_back(1);
  D4[2] = m_coordinates[2];
  D4[2].push_back(1);
  D4[3] = point;
  D4[3].push_back(1);

  double det0 = determinant(D0);
  double bc1 = determinant(D1)/det0;
  double bc2 = determinant(D2)/det0;
  double bc3 = determinant(D3)/det0;
  double bc4 = determinant(D4)/det0;

  double epsilon_ = 0.00001;

  /* assert convexity */
  assert(bc1 + bc2 + bc3 + bc4 < 1.0 + epsilon_ and 
         bc1 + bc2 + bc3 + bc4 > 1.0 - epsilon_ );

  if (bc1>=-epsilon_ and bc2>=-epsilon_ and bc3>=-epsilon_ and bc4>=-epsilon_) {
    return true; 
  } else {
    return false;
  }
}

vector<double> MasterElement::barycentric_coor(vector<double> point){
  if (m_dimension == 2){
    return barycentric_coor2D(point);
  } else {
    return barycentric_coor3D(point);
  }
}

vector<double> MasterElement::barycentric_coor2D(vector<double> point){

  vector<double> bcoord(3);

  vector<double> x1 = m_coordinates[0];
  vector<double> x2 = m_coordinates[1];
  vector<double> x3 = m_coordinates[2];

  vector<vector<double>> T(2);
  T[0].resize(2);
  T[1].resize(2);
  T[0][0] = x1[0] - x3[0];
  T[0][1] = x2[0] - x3[0];
  T[1][0] = x1[1] - x3[1];
  T[1][1] = x2[1] - x3[1];
  double detInv = determinant(T);
  bcoord[0] = ((x2[1] - x3[1])*(point[0] - x3[0]) + (x3[0] - x2[0])*(point[1] - x3[1]))/detInv;
  bcoord[1] = ((x3[1] - x1[1])*(point[0] - x3[0]) + (x1[0] - x3[0])*(point[1] - x3[1]))/detInv;
  bcoord[2] = 1.0 - bcoord[0] - bcoord[1];

  return bcoord;
}

vector<double> MasterElement::barycentric_coor3D(vector<double> point){

  vector<vector<double>> D0(4);
  D0[0] = m_coordinates[0];
  D0[0].push_back(1);
  D0[1] = m_coordinates[1];
  D0[1].push_back(1);
  D0[2] = m_coordinates[2];
  D0[2].push_back(1);
  D0[3] = m_coordinates[3];
  D0[3].push_back(1);

  vector<vector<double>> D1(4);
  D1[0] = point;
  D1[0].push_back(1);
  D1[1] = m_coordinates[1];
  D1[1].push_back(1);
  D1[2] = m_coordinates[2];
  D1[2].push_back(1);
  D1[3] = m_coordinates[3];
  D1[3].push_back(1);

  vector<vector<double>> D2(4);
  D2[0] = m_coordinates[0];
  D2[0].push_back(1);
  D2[1] = point;
  D2[1].push_back(1);
  D2[2] = m_coordinates[2];
  D2[2].push_back(1);
  D2[3] = m_coordinates[3];
  D2[3].push_back(1);

  vector<vector<double>> D3(4);
  D3[0] = m_coordinates[0];
  D3[0].push_back(1);
  D3[1] = m_coordinates[1];
  D3[1].push_back(1);
  D3[2] = point;
  D3[2].push_back(1);
  D3[3] = m_coordinates[3];
  D3[3].push_back(1);

  vector<vector<double>> D4(4);
  D4[0] = m_coordinates[0];
  D4[0].push_back(1);
  D4[1] = m_coordinates[1];
  D4[1].push_back(1);
  D4[2] = m_coordinates[2];
  D4[2].push_back(1);
  D4[3] = point;
  D4[3].push_back(1);
  
  double det0 = determinant(D0);
  double det1 = determinant(D1);
  double det2 = determinant(D2);
  double det3 = determinant(D3);
  double det4 = determinant(D4);

//  double epsilon_ = 0.000001;
//  assert(((det0 < det1 + det2 + det3 + det4 + epsilon_) or (det0 > det1 + det2 + det3 + det4 - epsilon_)) and
//          "Only acceptable to compute barycentric coordinates of a point inside the element.");
  assert(det0 != 0.0);
  vector<double> bcoor(4);

  bcoor[0] = det1/det0;
  bcoor[1] = det2/det0;
  bcoor[2] = det3/det0;
  bcoor[3] = det4/det0;

  return bcoor;
}

double MasterElement::a_grad_phi_i__dot__a_grad_phi_j(int i, int j, vector<double> a){
  return dot(transpose(m_J_inv) * grad_phi(i), a) * dot(transpose(m_J_inv) * grad_phi(j), a) * m_weights_sum * m_detJac;
}

//double MasterElement::a_grad_phi_i__dot__a_grad_phi_j(int i, int j, vector<double> a){
//  return dot(transpose(m_J_inv) * grad_phi(i), a) * dot(transpose(m_J_inv) * grad_phi(j), a) * m_weights_sum * m_detJac;
//}

//double MasterElement::a_grad_phi_i__dot__grad_phi_j(int i, int j, vector<double> a){
//  return dot(transpose(m_J_inv) * grad_phi(i), a) * transpose(m_J_inv) * grad_phi(j) * m_weights_sum * m_detJac;
//}

/* Transform quadrature rules from [-1,1] to [0, 1] (or equivalently for higher dimensions) */
//void MasterElement::arrangeReferenceDomain(){
//  if (m_dimension == 1){
//    for (int i = 0; i < m_quadraturePoints.size(); i++){
//      m_quadraturePoints[i][0] = 0.5 * m_quadraturePoints[i][0] + 0.5; 
//      m_weights[i] = 0.5 * m_weights[i];
//    }
//  } 
//}
