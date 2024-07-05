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

#include<masterElementBD.hpp>

void MasterElementBD::initialize(Parameters parameters, int dimension){

  par = parameters;
  m_dimension = dimension;
  m_nbNodesPerElement = m_dimension;
  m_nbDofsPerNode = par.nbDofsPerNode()[0];

  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "Element BD: Initializing." << endl;

  m_quadraturePoints.resize(par.nbQuadraturePointsBD());
  m_weights.resize(par.nbQuadraturePointsBD());
  for (int i = 0; i < par.nbQuadraturePointsBD(); i++){ 
    m_quadraturePoints[i].resize(m_dimension);  
  }

  if (m_dimension == 2){
    if (par.nbQuadraturePointsBD() == 1){

      m_quadraturePoints[0][0] = 0.5;
      m_quadraturePoints[0][1] = 0.0;
      m_quadraturePoints[0][2] = 0.0;
      m_weights[0] = 1.0;
    } else if (par.nbQuadraturePointsBD() == 3){
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
      errorMessage("MasterElementBD::initialize", "Quadrature rule not valid.");
    }
  }

  if (m_dimension == 3){
    if(par.nbQuadraturePointsBD() == 1){

      m_quadraturePoints[0][0] = 0.333333333333333333;
      m_quadraturePoints[0][1] = 0.333333333333333333;
      m_quadraturePoints[0][2] = 0.0;
      m_weights[0] = 0.5;

    } else if (par.nbQuadraturePointsBD() == 3){

      m_quadraturePoints[0][0] = 1.0/6.0;
      m_quadraturePoints[0][1] = 1.0/6.0;
      m_quadraturePoints[0][2] = 0.0;
      m_weights[0] = 1.0/6.0;

      m_quadraturePoints[1][0] = 1.0/6.0;
      m_quadraturePoints[1][1] = 2.0/3.0;
      m_quadraturePoints[1][2] = 0;
      m_weights[1] = 1.0/6.0;

      m_quadraturePoints[2][0] = 2.0/3.0;
      m_quadraturePoints[2][1] = 1.0/6.0;
      m_quadraturePoints[2][2] = 0.0;
      m_weights[2] = 1.0/6.0;

    } else if (par.nbQuadraturePointsBD() == 4){


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

    } else if(par.nbQuadraturePointsBD() == 6){
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
      errorMessage("MasterElementBD::initialize", "Quadrature rule not valid.");
    }
  }

  /* allocate elementary matrix and jacobian */
  m_J.resize(3);
  for (int i = 0; i < 3; i++){
    m_J[i].resize(m_dimension - 1);
  }
}

void MasterElementBD::setCoordinates(const vector<vector<double>> & coordinates){
  m_coordinates = coordinates;

  if (m_dimension == 2){
    m_J[0][0] = m_coordinates[1][0] - m_coordinates[0][0]; 
    m_J[1][0] = m_coordinates[1][1] - m_coordinates[0][1];
    m_J[2][0] = 0.0;
  }

  if (m_dimension == 3){
//    m_J[0][0] = m_coordinates[1][0] - m_coordinates[0][0]; 
//    m_J[1][0] = m_coordinates[1][1] - m_coordinates[0][1];
//    m_J[2][0] = m_coordinates[1][2] - m_coordinates[0][2];
//    m_J[0][1] = m_coordinates[2][0] - m_coordinates[0][0];
//    m_J[1][1] = m_coordinates[2][1] - m_coordinates[0][1];
//    m_J[2][1] = m_coordinates[2][2] - m_coordinates[0][2];
    m_J[0][0] = m_coordinates[1][0] - m_coordinates[0][0]; 
    m_J[1][0] = m_coordinates[1][1] - m_coordinates[0][1];
    m_J[2][0] = m_coordinates[1][2] - m_coordinates[0][2];

    m_J[0][1] = m_coordinates[2][0] - m_coordinates[0][0];
    m_J[1][1] = m_coordinates[2][1] - m_coordinates[0][1];
    m_J[2][1] = m_coordinates[2][2] - m_coordinates[0][2];
  }
  m_sqrtDetJTJ = sqrt( fabs(determinant(transpose(m_J)*m_J)) );
}

void MasterElementBD::computeVolume(){
  if (m_dimension == 3){
    vector<double> v1 = (m_coordinates[1] - m_coordinates[0]);
    vector<double> v2 = (m_coordinates[2] - m_coordinates[0]);
    m_volume = norm(cross(v1,v2))/2;
  } else if  (m_dimension == 2){
    m_volume = norm(m_coordinates[1] - m_coordinates[0]);
  }
}

void MasterElementBD::computeNormal(){
  m_normal.resize(3);
  if (m_dimension == 3){
    /* exterior normal if convention is well-followed on geometry file */
    vector<double> v1 = m_coordinates[1] - m_coordinates[0];
    vector<double> v2 = m_coordinates[2] - m_coordinates[0];
    m_normal = cross(v1,v2);
  } else if  (m_dimension == 2){
    m_normal[0] = -m_coordinates[0][1] + m_coordinates[1][1];
    m_normal[1] = -m_coordinates[1][0] + m_coordinates[0][0];
    m_normal[2] = 0.0;
  }
  m_normal = m_normal/norm(m_normal);
}

double MasterElementBD::mass(int i, int j){
  double mass = 0.0;
  for (int qp = 0; qp < par.nbQuadraturePointsBD(); qp++){
    mass = mass + phi(qp,i) * phi(qp,j) * m_weights[qp];
  }
  return mass * m_sqrtDetJTJ;
}

double MasterElementBD::phi(const int & qp_point_id, int id_shape_function){
  double phi_eval;
  if (m_dimension == 3){
    if (id_shape_function == 0){
      phi_eval = 1.0 - m_quadraturePoints[qp_point_id][0] - m_quadraturePoints[qp_point_id][1]; 
    } else if (id_shape_function == 1){
      phi_eval = m_quadraturePoints[qp_point_id][0];
    } else if (id_shape_function == 2){
      phi_eval = m_quadraturePoints[qp_point_id][1];
    }
  } else if (m_dimension == 2){
    if (id_shape_function == 0){
      phi_eval = 1.0 - m_quadraturePoints[qp_point_id][0]; 
    } else if (id_shape_function == 1){
      phi_eval = m_quadraturePoints[qp_point_id][0];
    }
  }
  return phi_eval; 
}

double MasterElementBD::computeFlow(vector<double> vectorField){
  double integral = 0.0;
  for (int i = 0; i < m_nbNodesPerElement; i++){
    for (int comp = 0; comp < m_nbDofsPerNode; comp++){
      for (int qp = 0; qp < par.nbQuadraturePointsBD(); qp++){
        integral += vectorField[m_nbDofsPerNode*i+comp] * m_normal[comp] * phi(qp, i) * m_weights[qp] * fabs(m_sqrtDetJTJ);
      }
    }
  }
  return integral;
}

//void MasterElementBD::arrangeReferenceDomain(){
//  if (m_dimension == 1){
//    for (int i = 0; i < m_quadraturePoints.size(); i++){
//      m_quadraturePoints[i][0] = 0.5 * m_quadraturePoints[i][0] + 0.5; 
//      m_weights[i] = 0.5 * m_weights[i];
//    }
//  } 
//}
