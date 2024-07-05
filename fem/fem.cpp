/*=============================================================================
  This file is part of the code MAD
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

#include<fem.hpp>

void FEM::initialize(Parameters parameters, const Geometry & geometry, const Boundary & boundary, Mat & matrix){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "FEM: Initializing" << endl;

  par = parameters;
  geo = geometry;
  bd = boundary;

  fe.initialize(par, geo.dimension());
  feBD.initialize(par, geo.dimension());
  m_verbose = par.verbose();

  nbDofVar.resize(par.nbVariables());

  nbDofs = 0;
  for (int i = 0; i < par.nbDofsPerNode().size(); i++){
    nbDofs += par.nbDofsPerNode()[i]*geo.nbVertices;
    nbDofVar[i] = par.nbDofsPerNode()[i]*geo.nbVertices;
  }
  nbVertices = geo.nbVertices;
  nbNodesPerElement = geo.dimension()+1; 

  A = matrix;
  code = MatGetOwnershipRange(A, &m, &n); CHKERR(code);

}

void FEM::a_grad_u_dot_v(double coef, vector<double> u_el, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  for (int j = 0; j < nbNodesPerElement; j++){
    vector<double> u_node(geo.dimension());
    for (int comp = 0; comp < geo.dimension(); comp++){
      u_node[comp] = u_el[geo.dimension()*j+comp];
    }
    for (int i = 0; i < nbNodesPerElement; i++){
      for (int k = 0; k < par.nbDofsPerNode()[testLabel]; k++){
        int row = par.nbDofsPerNode()[testLabel]*m_simplex[i]+k;
        vector<double> conv = fe.dphi_i_phi_j_phi_i(j,i);
//        vector<double> conv = fe.dphi_i_phi_j(j,i);
        if (row >= m && row < n){
          for (int l = 0; l < par.nbDofsPerNode()[trialLabel]; l++){
            int col = par.nbDofsPerNode()[trialLabel]*m_simplex[j]+l;
            code = MatSetValue(A, row, col, u_node[l]*coef*conv[l], ADD_VALUES); CHKERR(code);
          }
        }
      }
    }
  }
}

void FEM::u__dot__a_grad_v(double coef, vector<double> u_el, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  for (int j = 0; j < nbNodesPerElement; j++){
    vector<double> u_node(geo.dimension());
    for (int comp = 0; comp < geo.dimension(); comp++){
      u_node[comp] = u_el[geo.dimension()*j+comp];
    }
    for (int i = 0; i < nbNodesPerElement; i++){
      for (int k = 0; k < par.nbDofsPerNode()[testLabel]; k++){
        int row = par.nbDofsPerNode()[testLabel]*m_simplex[i]+k;
        vector<double> conv = fe.dphi_i_phi_j(i,j);
        if (row >= m && row < n){
          for (int l = 0; l < par.nbDofsPerNode()[trialLabel]; l++){
            int col = par.nbDofsPerNode()[trialLabel]*m_simplex[j]+l;
            code = MatSetValue(A, row, col, u_node[l]*coef*conv[l], ADD_VALUES); CHKERR(code);
          }
        }
      }
    }
  }
}


void FEM::a_grad_u__dot__a_grad_v(double coef, vector<double> u_el, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  for (int j = 0; j < nbNodesPerElement; j++){
    vector<double> u_node(geo.dimension());
    for (int comp = 0; comp < geo.dimension(); comp++){
      u_node[comp] = u_el[geo.dimension()*j+comp];
    }
    for (int i = 0; i < nbNodesPerElement; i++){
      double Aij = fe.stiff(i,j);
      for (int k = 0; k < par.nbDofsPerNode()[testLabel]; k++){
        int row = par.nbDofsPerNode()[testLabel]*m_simplex[i]+k;
        if (row >= m && row < n){
          for (int l = 0; l < par.nbDofsPerNode()[testLabel]; l++){
            int col = par.nbDofsPerNode()[trialLabel]*m_simplex[j]+l;
            code = MatSetValue(A, row, col , coef*u_node[k]*u_node[l]*Aij, ADD_VALUES); CHKERR(code);
          }
        }
      }
    }
  }
}

void FEM::a_grad_u__dot__grad_q(double coef, vector<double> u_el, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  for (int j = 0; j < nbNodesPerElement; j++){
    vector<double> u_node(geo.dimension());
    for (int comp = 0; comp < geo.dimension(); comp++){
      u_node[comp] = u_el[geo.dimension()*j+comp];
    }
    for (int i = 0; i < nbNodesPerElement; i++){
      for (int comp = 0; comp < par.nbDofsPerNode()[trialLabel]; comp++){
        double Aij = coef * fe.stiff(i,j) * u_node[comp];
        int row = m_simplex[i];
        if (row >= m && row < n){
          int col = par.nbDofsPerNode()[trialLabel]*m_simplex[j]+comp;
          code = MatSetValue(A, row, col , Aij, ADD_VALUES); CHKERR(code);
        }
      }
    }
  }
}

void FEM::a_grad_v__dot__grad_p(double coef, vector<double> u_el, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  for (int j = 0; j < nbNodesPerElement; j++){
    vector<double> u_node(geo.dimension());
    for (int comp = 0; comp < geo.dimension(); comp++){
      u_node[comp] = u_el[geo.dimension()*j+comp];
    }
    for (int i = 0; i < nbNodesPerElement; i++){
      for (int comp = 0; comp < par.nbDofsPerNode()[testLabel]; comp++){
        double Aij = coef * fe.stiff(i,j) * u_node[comp];
        int row = par.nbDofsPerNode()[testLabel]*m_simplex[i]+comp;
        if (row >= m && row < n){
          int col = m_simplex[j];
          code = MatSetValue(A, row, col , Aij, ADD_VALUES); CHKERR(code);
        }
      }
    }
  }
}


void FEM::epsilon_u_epsilon_v(double coef, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  /* Assemble elementary matrices */
  for (int j = 0; j < nbNodesPerElement; j++){
    for (int i = 0; i < nbNodesPerElement; i++){

      /* diagonal components */
      for (int comp = 0; comp < par.nbDofsPerNode()[testLabel]; comp++){
        int row = offset_row + par.nbDofsPerNode()[testLabel]*m_simplex[i]+comp;
        int col = offset_col + par.nbDofsPerNode()[trialLabel]*m_simplex[j]+comp;
        if (row >= m && row < n){
          double Kij = fe.stiff(i,j); 
          code = MatSetValue(A, row, col, coef*Kij, ADD_VALUES); CHKERR(code);
        }
      }

      /* deviatoric components */
      for (int comp_v = 0; comp_v < par.nbDofsPerNode()[testLabel]; comp_v++){
        int row = offset_row + par.nbDofsPerNode()[testLabel]*m_simplex[i]+comp_v;
        if (row >= m && row < n){
          for (int comp_u = 0; comp_u < par.nbDofsPerNode()[trialLabel]; comp_u++){
            int col = offset_col + par.nbDofsPerNode()[trialLabel]*m_simplex[j]+comp_u;
            double Kdev = 0.5*fe.dphi_dx(i, comp_u) * fe.dphi_dx(j, comp_v) * fe.detJacobian() * fe.weightsSum();
            code = MatSetValue(A, row, col,  coef*Kdev, ADD_VALUES); CHKERR(code);
          }
        }
      }
    }
  }
}

void FEM::du_x_dv_y(double coef, int x, int y, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  /* Assemble elementary matrices */
  for (int i = 0; i < nbNodesPerElement; i++){
    for (int j = 0; j < nbNodesPerElement; j++){
      double Kij = fe.dphi_i_COMP_dphi_j_COMP(i,j,x,y); 
      int row = offset_row + m_simplex[i];
      int col = offset_col + m_simplex[j];
      if (row >= m && row < n){
        code = MatSetValue(A, row, col, coef*Kij, ADD_VALUES); CHKERR(code);
      }
    }
  }
}

void FEM::u_dv_y(double coef, int y, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  /* Assemble elementary matrices */
  for (int i = 0; i < nbNodesPerElement; i++){
    for (int j = 0; j < nbNodesPerElement; j++){
      double Kij = fe.phi_i_dphi_j_COMP(j,i,y); 
      int row = offset_row + m_simplex[i];
      int col = offset_col + m_simplex[j];
      if (row >= m && row < n){
        code = MatSetValue(A, row, col, coef*Kij, ADD_VALUES); CHKERR(code);
      }
    }
  }
}

void FEM::a_du_x_v(double coef, int x, vector<double> a, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  /* Assemble elementary matrices */
  for (int i = 0; i < nbNodesPerElement; i++){
    for (int j = 0; j < nbNodesPerElement; j++){
      double Kij = fe.phi_dphi_i_COMP_phi_j(j,i,x,a); 
      int row = offset_row + m_simplex[i];
      int col = offset_col + m_simplex[j];
      if (row >= m && row < n){
        code = MatSetValue(A, row, col, coef*Kij, ADD_VALUES); CHKERR(code);
      }
    }
  }
}

void FEM::a_du_x_b_dv_y(double coef, int x, int y, vector<double> a, vector<double> b, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  /* Assemble elementary matrices */
  for (int i = 0; i < nbNodesPerElement; i++){
    for (int j = 0; j < nbNodesPerElement; j++){
      double Kij = fe.phi_phi_dphi_i_COMP_dphi_j_COMP(i,j,x,y,a,b); 
      int row = offset_row + m_simplex[i];
      int col = offset_col + m_simplex[j];
      if (row >= m && row < n){
        code = MatSetValue(A, row, col, coef*Kij, ADD_VALUES); CHKERR(code);
      } /* this exploits the symmetry of the bi-linear form */
      if (col >= m && col < n){
        code = MatSetValue(A, col, row, coef*Kij, ADD_VALUES); CHKERR(code);
      }
    }
  }
}

void FEM::du_x_v(double coef, int x, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  /* Assemble elementary matrices */
  for (int i = 0; i < nbNodesPerElement; i++){
    for (int j = 0; j < nbNodesPerElement; j++){
      double Kij = fe.phi_i_dphi_j_COMP(i,j,x); 
      int row = offset_row + m_simplex[i];
      int col = offset_col + m_simplex[j];
      if (row >= m && row < n){
        code = MatSetValue(A, row, col, coef*Kij, ADD_VALUES); CHKERR(code);
      }
    }
  }
}

void FEM::u_dot_v_scalar(double coef, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  /* Assemble elementary matrices */
  for (int i = 0; i < nbNodesPerElement; i++){
    for (int j = 0; j < nbNodesPerElement; j++){
      double Mij = fe.phi_i_phi_j(i,j);
      int row = offset_row + m_simplex[i];
      int col = offset_col + m_simplex[j];
      if (row >= m && row < n){
        code = MatSetValue(A, row, col, coef*Mij, ADD_VALUES); CHKERR(code);
      }
    }
  }
}

void FEM::u_dot_v(double coef, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  /* Assemble elementary matrices */
  for (int i = 0; i < nbNodesPerElement; i++){
    for (int j = 0; j < nbNodesPerElement; j++){
      double Mij = fe.mass(i,j);
      for (int comp = 0; comp < par.nbDofsPerNode()[testLabel]; comp++){
        int row = offset_row + par.nbDofsPerNode()[testLabel]*m_simplex[i]+comp;
        int col = offset_col + par.nbDofsPerNode()[trialLabel]*m_simplex[j]+comp;
        if (row >= m && row < n){
          code = MatSetValue(A, row, col, coef*Mij, ADD_VALUES); CHKERR(code);
        }
      }
    }
  }
}

void FEM::grad_u_grad_v(double coef, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  /* Assemble elementary matrices */
  for (int i = 0; i < nbNodesPerElement; i++){
    for (int j = 0; j < nbNodesPerElement; j++){
      double Kij = fe.stiff(i,j); 
      for (int comp = 0; comp < par.nbDofsPerNode()[testLabel]; comp++){
        int row = offset_row + par.nbDofsPerNode()[testLabel]*m_simplex[i]+comp;
        int col = offset_col + par.nbDofsPerNode()[trialLabel]*m_simplex[j]+comp;
        if (row >= m && row < n){
          code = MatSetValue(A, row, col, coef*Kij, ADD_VALUES); CHKERR(code);
        }
      }
    }
  }
}

void FEM::div_u_v(double coef, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  /* Assemble elementary matrices */
  for (int j = 0; j < nbNodesPerElement; j++){
    for (int i = 0; i < nbNodesPerElement; i++){
      for (int comp = 0; comp < par.nbDofsPerNode()[trialLabel]; comp++){
        int row = offset_row + m_simplex[i];
        int col = offset_col + par.nbDofsPerNode()[trialLabel]*m_simplex[j]+comp;
        if (row  >= m && row < n){
          double Bij = fe.mixed(j,i,comp);
          code = MatSetValue(A, row, col, coef*Bij, ADD_VALUES); CHKERR(code); /* B^T */ 
        }
      }
    }
  }
}

void FEM::div_u_div_v(double coef, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  for (int j = 0; j < nbNodesPerElement; j++){
    for (int i = 0; i < nbNodesPerElement; i++){
      for (int comp_v = 0; comp_v < par.nbDofsPerNode()[testLabel]; comp_v++){
        /* deviatoric components */
        int row = par.nbDofsPerNode()[testLabel]*m_simplex[i]+comp_v;
        if (row >= m && row < n){
          for (int comp_u = 0; comp_u < par.nbDofsPerNode()[trialLabel]; comp_u++){
            double div_div_ij = fe.dphi_dx(i, comp_v) * fe.dphi_dx(j, comp_u) * fe.detJacobian()*fe.weightsSum();
            int col = par.nbDofsPerNode()[trialLabel]*m_simplex[j]+comp_u;
            code = MatSetValue(A, row, col, div_div_ij, ADD_VALUES); CHKERR(code);
          }
        }
      }
    }
  }
}

void FEM::u_div_v(double coef, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  /* Assemble elementary matrices */
  for (int j = 0; j < nbNodesPerElement; j++){
    for (int i = 0; i < nbNodesPerElement; i++){
      for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
        int row = offset_row + par.nbDofsPerNode()[testLabel]*m_simplex[i]+comp;
        int col = offset_col + m_simplex[j];
        if (row >= m && row < n){
          double Bij = fe.mixed(i,j,comp);
          code = MatSetValue(A, row, col, coef*Bij, ADD_VALUES); CHKERR(code); /* -B*/ 
        }
      }
    }
  }
}

void FEM::u__dot__grad_q(double coef, int testLabel, int trialLabel){

  int offset_row = 0;
  for (int i = 0; i < testLabel; i++){
    offset_row += nbDofVar[i];
  }

  int offset_col = 0;
  for (int i = 0; i < trialLabel; i++){
    offset_col += nbDofVar[i];
  }

  for (int j = 0; j < nbNodesPerElement; j++){
    for (int i = 0; i < nbNodesPerElement; i++){
      for (int comp = 0; comp < par.nbDofsPerNode()[trialLabel]; comp++){
        vector<double> Aij = coef * fe.dphi_i_phi_j(i,j);
        int row = m_simplex[i];
        if (row >= m && row < n){
          int col = par.nbDofsPerNode()[trialLabel]*m_simplex[j]+comp;
          code = MatSetValue(A, row, col , Aij[comp], ADD_VALUES); CHKERR(code);
        }
      }
    }
  }
}

void FEM::setSimplex(vector<int> element){
  m_simplex = element;
  /* get finite element coordinates */
  vector<vector<double>> coordinates(nbNodesPerElement);
  for (int i = 0; i < nbNodesPerElement; i++){
    coordinates[i] = geo.coordinates()[m_simplex[i]];
  }
  fe.setCoordinates(coordinates);
  fe.computeSize();
}

void FEM::copySimplex(FEM & finiteElement){
  m_simplex = finiteElement.simplex(); 
  fe = finiteElement.finiteElement();
}
