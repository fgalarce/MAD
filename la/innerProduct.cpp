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

#include<innerProduct.hpp>

void InnerProduct::finalize(){
  MatDestroy(&matrix);
  if (par.innerProduct() == "H1"){
    MatDestroy(&m_stiff);
  }
}

///* old initializer */
//void InnerProduct::initialize(Parameters parameters, Geometry geo, const Mat & M){
//  cout << "IP: Initializing inner product." << endl;
//  par = parameters;
//  m_matSize = geo.nbVertices*geo.dimension();
//  m_nbVertices = geo.nbVertices;
//  if (par.mixed_problem()){
//    m_matSize = m_matSize + geo.nbVertices;
//  }
//  MatDuplicate(M, MAT_COPY_VALUES, &matrix);
//}

///* old initializer */
//void InnerProduct::initialize(Parameters parameters, const Geometry & geo){
//
//  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
//  if (m_world_rank == 0) cout << "IP: Initializing inner product." << endl;
//  par = parameters;
//
//  assert(par.innerProduct() == "L2" || par.innerProduct() == "l2" || par.innerProduct() == "H1");
//  if (par.nbVariables() > 1){
//    assert( (par.innerProduct(0) == "H1" && par.innerProduct(1) == "L2")
//         || (par.innerProduct(0) == "L2" && par.innerProduct(1) == "L2")
//         || (par.innerProduct(0) == "l2" && par.innerProduct(1) == "l2"));
//  }
//
//  /* Ambient space matrix */
//  m_nbDofsPerNode.resize(par.nbVariables());
//  m_matSize = 0;
//  for (int i = 0; i < par.nbVariables(); i++){
//    m_matSize = m_matSize + geo.nbVertices*par.nbDofsPerNode()[i];
//    m_nbDofsPerNode[i] = par.nbDofsPerNode()[i];
//  }
//  m_nbVertices = geo.nbVertices;
//  m_geo = geo;
//  mat(matrix, m_matSize, m_matSize);
//  mat(m_stiff, m_matSize, m_matSize);
//
//  /* build inner product matrix if not provided in parameters file */
//  if (par.mass_matrix_path() == ""  && par.innerProduct() != "l2"){ 
//    assembleAmbientSpaceMatrix();
//
//  /* otherwise load them from file */
//  } else { 
//    loadAmbientSpaceMatrix();
//  }
//
//  if (par.innerProduct() != "l2"){
//    PetscBool isSymmetric;
//    code = MatIsSymmetric(matrix, 0.00000001, &isSymmetric); CHKERR(code);
//    assert(isSymmetric && "IP: matrix is not symmetric.");
//  }
//}

void InnerProduct::initialize(Parameters parameters, const Geometry & geometry, const Boundary & boundary){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  par = parameters;
  geo = geometry;
  bd = boundary;

  nbDofVar.resize(par.nbVariables());

  nbDofs = 0;
  for (int i = 0; i < par.nbDofsPerNode().size(); i++){
    nbDofs += par.nbDofsPerNode()[i]*geo.nbVertices;
    nbDofVar[i] = par.nbDofsPerNode()[i]*geo.nbVertices;
  }
  nbVertices = geo.nbVertices;
  nbNodesPerElement = geo.dimension()+1; 

  fe.initialize(par, geo.dimension());
  mat(m_ip, nbDofs, nbDofs);
  code = MatGetOwnershipRange(m_ip, &m, &n); CHKERR(code);

  fem.initialize(par, geo, bd, m_ip);
  assembleInnerProduct();
}

double InnerProduct::operator()(const Vec & u, const Vec & v){
  double ip;
  Vec aux = vec(nbDofs);
  code = MatMult(m_ip, u, aux); CHKERR(code);
  code = VecDot(aux, v, &ip); CHKERR(code);
  VecDestroy(&aux);
  return ip;
}

void InnerProduct::assembleInnerProduct(){
  if (m_world_rank == 0) cout << "IP: Assembling inner product matrix." << endl;
  for (int partId = 0; partId < geo.elements().size(); partId++){
    for (vector<int> simplex : geo.elements()[partId]){ /* loop on elements */
      fem.setSimplex(simplex);
      for (int idVar = 0; idVar < par.nbVariables(); idVar++){
        if (par.innerProduct(idVar) == "L2" or par.innerProduct(idVar) == "H1"){
          fem.u_dot_v(par.factorInnerProduct()[idVar], idVar, idVar);
        }
        if (par.innerProduct(idVar) == "H1"){
            fem.grad_u_grad_v(par.factorInnerProduct()[idVar], idVar, idVar);
        }
      } 
    }
  }
  /* fill with identities wherever l2 metric is chosen */
  for (int i = 0; i < par.nbVariables(); i++){
    if (par.innerProduct(i) == "l2"){
      for (int j = i*geo.nbVertices; j < (i+1)*geo.nbVertices; j++){
        matSet(m_ip, j, j, 1.0);
      }
    }
  }
  code = MatAssemblyBegin(m_ip, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(m_ip, MAT_FINAL_ASSEMBLY); CHKERR(code);
}

double InnerProduct::operator()(const Vec & u, const Vec & v, int varLabel){

  double ip;
  int vecSize, vecSize_v;
  VecGetSize(u, &vecSize);
  VecGetSize(v, &vecSize_v);
  assert(vecSize == vecSize_v);

  if (par.innerProduct() == "l2"){
    code = VecDot(u, v, &ip); CHKERRQ(code);

  } else if (par.innerProduct() != "l2" && varLabel == -1){
    Vec mtv;
    vec(mtv, m_matSize);
    code = MatMult(matrix, v, mtv); CHKERRQ(code);
    code = VecDot(u, mtv, &ip); CHKERRQ(code);
    code = VecDestroy(&mtv); CHKERRQ(code);

  } else if (par.innerProduct() != "l2" && varLabel != -1) {

    Vec u_mixed = zeros(m_matSize);  
    Vec v_mixed = zeros(m_matSize);  

    int offset = 0;
    for (int i = 0; i < varLabel; i++){
      offset = offset + par.nbDofsPerNode()[i] * m_nbVertices;
    }
    
    int nbDofs = par.nbDofsPerNode()[varLabel]*m_nbVertices;
    int m,n;
    code = VecGetOwnershipRange(u, &m, &n); CHKERR(code);
    for (int i = offset; i < offset + nbDofs; i++){
      if (i >= m && i < n){
        double var_value;
        code = VecGetValues(u, 1, &i, &var_value); CHKERRQ(code); 
        code = VecSetValue(u_mixed, i, var_value, INSERT_VALUES); CHKERRQ(code); 
        code = VecGetValues(v, 1, &i, &var_value); CHKERRQ(code); 
        code = VecSetValue(v_mixed, i, var_value, INSERT_VALUES); CHKERRQ(code); 
      }
    }
    code = VecAssemblyBegin(u_mixed); CHKERRQ(code);
    code = VecAssemblyEnd(u_mixed); CHKERRQ(code);
    code = VecAssemblyBegin(v_mixed); CHKERRQ(code);
    code = VecAssemblyEnd(v_mixed); CHKERRQ(code);

    Vec mtv;
    vec(mtv, m_matSize);
    code = MatMult(matrix, v_mixed, mtv); CHKERRQ(code);
    code = VecDot(u_mixed, mtv, &ip); CHKERRQ(code);
    code = VecDestroy(&mtv); CHKERRQ(code);

  } else {
    exit(1);
    cout << "WRONG" << endl;
  }

  return ip;
}

double InnerProduct::operator()(const Vec & u, const Mat & M, const Vec & v){
  int matSize;
  code = MatGetSize(M, &matSize, NULL); CHKERR(code);
  double ip;
  Vec mtv;
  vec(mtv, matSize);
  code = MatMult(M, v, mtv); CHKERRQ(code);
  code = VecDot(u, mtv, &ip); CHKERRQ(code);
  code = VecDestroy(&mtv); CHKERRQ(code);
  return ip;
}

double InnerProduct::only_grad(const Vec & u, const Vec & v){
  PetscReal ip;
  Vec mtv;
  vec(mtv, m_matSize);
  code = MatMult(m_stiff, v, mtv); CHKERR(code);
  code = VecDot(u, mtv, &ip); CHKERR(code);
  code = VecDestroy(&mtv); CHKERR(code);
  return ip;
}

/* old function kept for compatibility */
void InnerProduct::assembleAmbientSpaceMatrix(){

  if (m_world_rank == 0) cout << "IP: Assembling ambient space matrix\n";
  int nbNodesPerElement = m_geo.dimension()+1; 

  MasterElement fe; 
  fe.initialize(par, m_geo.dimension());

  int m,n;
  code = MatGetOwnershipRange(matrix, &m, &n); CHKERR(code);
  int nbDofs0 = m_nbDofsPerNode[0]*m_nbVertices;

  int feId = 0;
  for (int partId = 0; partId < m_geo.elements().size(); partId++){ /* loop on parts */

    for (vector<int> simplex : m_geo.elements()[partId]){ /* loop on elements */
      vector<int> elementNodes = m_geo.elements()[partId][feId]; feId++;
      vector<vector<double>> coordinates(nbNodesPerElement);
      /* get finite element coordinates */
      for (int nodeId = 0; nodeId < nbNodesPerElement; nodeId++){
        coordinates[nodeId] = m_geo.coordinates()[elementNodes[nodeId]];
      }
      fe.setCoordinates(coordinates);

      /* set values from local to global */
      for (int i = 0; i < nbNodesPerElement; i++){
        for (int j = 0; j < nbNodesPerElement; j++){
          for (int comp = 0; comp < m_nbDofsPerNode[0]; comp++){
            if ((m_nbDofsPerNode[0]*simplex[i]+comp) >= m && (m_nbDofsPerNode[0]*simplex[i]+comp) < n){
              double Mij = fe.mass(i,j);
              double Kij = fe.stiffness(i,j);
              double Matrix_ij;
              if (par.innerProduct(0) == "H1"){
                Matrix_ij = Mij + par.weight_stiffness() * Kij;
              } else {
                Matrix_ij = Mij;
              }
              code = MatSetValue(matrix, m_nbDofsPerNode[0]*simplex[i]+comp, m_nbDofsPerNode[0]*simplex[j]+comp, Matrix_ij, ADD_VALUES); CHKERR(code);
              code = MatSetValue(m_stiff, m_nbDofsPerNode[0]*simplex[i]+comp, m_nbDofsPerNode[0]*simplex[j]+comp, Kij, ADD_VALUES); CHKERR(code);
            }
          }
          /* mixed problems */
          if (par.nbVariables() > 1){
            for (int comp = 0; comp < m_nbDofsPerNode[1]; comp++){
              if ((m_nbDofsPerNode[1]*simplex[i]+comp) >= m && (m_nbDofsPerNode[1]*simplex[i]+comp) < n){
                /* mass */
                double Mij = fe.mass(i,j);
                /* stiffness */
                double Kij = fe.stiffness(i,j);
                double Matrix_ij;
                if (par.innerProduct(1) == "H1"){
                  Matrix_ij = par.weight_mass() * Mij + par.weight_stiffness() * Kij;
                } else {
                  Matrix_ij = par.weight_mass() * Mij;
                }
                code = MatSetValue(matrix, nbDofs0 + m_nbDofsPerNode[1]*simplex[i]+comp, nbDofs0 + m_nbDofsPerNode[1]*simplex[j]+comp, Matrix_ij, ADD_VALUES); CHKERR(code);
              }
            }
          }
        }
      }
    }
  }
  code = MatAssemblyBegin(m_stiff, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(m_stiff, MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY); CHKERR(code);
}

void InnerProduct::loadAmbientSpaceMatrix(){
  /* Import mass and stiffness matrices */
  if (par.innerProduct() == "L2"){
    loadMat(matrix, par.mass_matrix_path());
    if (m_world_rank == 0) cout << "IP: Frobenius norm of mass matrix: " << norm(matrix) << "." << endl;

  } else if (par.innerProduct() == "H1"){
    mat(m_stiff, m_matSize, m_matSize);
    loadMat(matrix, par.mass_matrix_path());
    loadMat(m_stiff, par.stiffness_matrix_path());
    code = MatAssemblyBegin(matrix, MAT_FINAL_ASSEMBLY); CHKERR(code); 
    code = MatAssemblyEnd(matrix, MAT_FINAL_ASSEMBLY); CHKERR(code);
    code = MatAssemblyBegin(m_stiff, MAT_FINAL_ASSEMBLY); CHKERR(code); 
    code = MatAssemblyEnd(m_stiff, MAT_FINAL_ASSEMBLY); CHKERR(code);

    code = MatAXPY(matrix, par.weight_stiffness(), m_stiff, SAME_NONZERO_PATTERN); /* M + K */
    if (m_world_rank == 0) cout << "IP: Frobenius norm of mass matrix: " << norm(matrix) << "." << endl;
    if (m_world_rank == 0) cout << "IP: Frobenius norm of stiffness matrix: " << norm(m_stiff) << "." << endl;
    if (m_world_rank == 0) cout << "IP: Frobenius norm of ambient space matrix: " << norm(matrix) << "." << endl;

  } else if (par.innerProduct() == "H1xL2"){
    m_matSize = m_matSize + m_geo.nbVertices;
    loadMat(matrix, par.mass_matrix_path());
    mat(m_stiff, m_matSize, m_matSize);
    loadMat(m_stiff, par.stiffness_matrix_path());
    code = MatAXPY(matrix, par.weight_stiffness(), m_stiff, SAME_NONZERO_PATTERN); /* M + K */
    if (m_world_rank == 0) cout << "IP: Frobenius norm of mass matrix: " << norm(matrix) << "." << endl;
    if (m_world_rank == 0) cout << "IP: Frobenius norm of stiffness matrix: " << norm(m_stiff) << "." << endl;
    if (m_world_rank == 0) cout << "IP: Frobenius of ambient space matrix: " << norm(matrix) << "." << endl;

    Mat mass_pressure;
    mat(mass_pressure, m_matSize, m_matSize);
    loadMat(mass_pressure, par.mass_matrix_path_pressure());
    MatAXPY(matrix, par.weight_pressure(), mass_pressure, SAME_NONZERO_PATTERN);
  }
}
