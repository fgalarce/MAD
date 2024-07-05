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

#include<linearAlgebra.hpp>

/* If not InnerProduct is passed, l2 space is assumed */
void LinearAlgebra:: initialize(Parameters parameters){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "LA: initializing l2." << endl;
  par = parameters;
}

void LinearAlgebra:: initialize(Parameters parameters, InnerProduct iproduct){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "LA: initializing." << endl;
  par = parameters;
  ip = iproduct;
}

bool LinearAlgebra::checkOrthonormality(vector<Vec> u){
  double epsilon_ = 0.00001;
  for (int i=0; i<u.size(); i++){
    for (int j=0; j<=i; j++){
      double uij = ip(u[i], u[j]);
      if ((uij <= delta(i,j) - epsilon_) or (uij >= delta(i,j) + epsilon_)){
        return false;
      }
    }
  }
  return true;
}

vector<double> LinearAlgebra::orthonormalize(vector<Vec> & u){
  if (m_world_rank == 0) cout << "LA: Orthonormalizing vectors with modified GS." << endl;
  vector<double> norms(u.size());
  for (int i = 0; i < u.size(); i++){
    loop_cout(i, u.size(), to_string(int(float(i)/float(u.size())*100)) + " %");
    for (int j = 0; j < i; j++){
      code = VecAXPY(u[i], -ip(u[i], u[j]), u[j]); CHKERR(code);
    }
    norms[i] = sqrt(ip(u[i], u[i]));
    VecScale(u[i], 1.0/norms[i]);
  }
  return norms;
}

vector<vector<double>> LinearAlgebra::eigenValueProblem(Mat & A, string mode){
  PetscErrorCode code;

  EPS eps;
  int nconv;
  EPSCreate(PETSC_COMM_WORLD, &eps);
  code = EPSSetOperators(eps, A, NULL); CHKERR(code);
  if (mode == "HEP"){
    code = EPSSetProblemType(eps, EPS_HEP); CHKERR(code);
  } else if (mode == "NHEP"){
    code = EPSSetProblemType(eps, EPS_GNHEP); CHKERR(code);
  } else {
    if (m_world_rank == 0) cout << "EPS: Eigen value problem type '" << mode << "' not valid.\n";
    exit(1);
  }

  code = EPSSetType(eps, EPSKRYLOVSCHUR); CHKERR(code);
  int nRows, nCols, dim;
  code = MatGetSize(A, &nRows ,&nCols); CHKERR(code);
  if (nRows > nCols){
    dim = nRows;  
  } else {
    dim = nCols;
  }
  code = EPSSetDimensions(eps, dim, PETSC_DECIDE, PETSC_DECIDE); CHKERR(code);
  code = EPSSetFromOptions(eps); CHKERR(code);
  code = EPSSolve(eps); CHKERR(code);
  code = EPSGetConverged(eps, &nconv); CHKERR(code);
  if (m_world_rank == 0) cout << "EPS: nconv " << nconv << endl; 

  vector<vector<double>> eValues(2); // a + bi
  eValues[0].resize(nconv);
  eValues[1].resize(nconv);
  for (int i=0; i<nconv; i++) {
    code = EPSGetEigenvalue(eps, i, &eValues[0][i], &eValues[1][i]); CHKERR(code);
  }

  EPSDestroy(&eps);
  return eValues;
}

vector<vector<double>> LinearAlgebra::eigenValueProblem(Mat & A, Mat & B, string mode){
  PetscErrorCode code;

  EPS eps;
  int nconv;
  EPSCreate(PETSC_COMM_WORLD, &eps);
  code = EPSSetOperators(eps, A, B); CHKERR(code);
  if (mode == "GNHEP"){
    code = EPSSetProblemType(eps, EPS_GNHEP); CHKERR(code);
  } else if (mode == "GHEP"){
    code = EPSSetProblemType(eps, EPS_GHEP); CHKERR(code);
  } else if (mode == "HEP"){
    code = EPSSetProblemType(eps, EPS_HEP); CHKERR(code);
  } else if (mode == "PGNHEP"){
    code = EPSSetProblemType(eps, EPS_PGNHEP); CHKERR(code);
  } else {
    if (m_world_rank == 0) cout << "EPS: Eigen value problem type '" << mode << "' not valid.\n";
    exit(1);
  }
  code = EPSSetType(eps, EPSKRYLOVSCHUR); CHKERR(code);
  int nRows, nCols, dim;
  code = MatGetSize(A, &nRows ,&nCols); CHKERR(code);
  if (nRows > nCols){
    dim = nRows;  
  } else {
    dim = nCols;
  }
  code = EPSSetDimensions(eps, dim, PETSC_DECIDE, PETSC_DECIDE); CHKERR(code);
  code = EPSSetFromOptions(eps); CHKERR(code);
  code = EPSSolve(eps); CHKERR(code);
  code = EPSGetConverged(eps, &nconv); CHKERR(code);
  if (m_world_rank == 0) cout << "EPS: nconv " << nconv << endl; 

  vector<vector<double>> eValues(2); // a + bi
  eValues[0].resize(nconv);
  eValues[1].resize(nconv);
  for (int i=0; i<nconv; i++) {
    code = EPSGetEigenvalue(eps, i, &eValues[0][i], &eValues[1][i]); CHKERR(code);
  }

  EPSDestroy(&eps);
  return eValues;
}

/* Returns (PTP)^{-1}PT or a good approximation of it */
Mat & LinearAlgebra::pseudoinverse(Mat A){

  PetscErrorCode code; 
  int m, n;
  MatGetSize(A, &m, &n);
 
  /* Singular value decomposition */ 
  int l = m;
  if (n < m){
    l = n;
  }
  Mat U = mat(m,l, "dense");
  Mat S = mat(l,l, "dense");
  Mat V = mat(n,l, "dense");
  SVD svd;           
  SVDCreate(PETSC_COMM_WORLD, &svd);
  SVDSetOperator(svd, A);
  SVDSetDimensions(svd, l, PETSC_DEFAULT, PETSC_DEFAULT);  
  SVDSetFromOptions(svd);
  SVDSolve(svd);
  for (int i = 0; i < l; i++){
    double svalue;
    double v_array[m];
    Vec v = vec(n);
    Vec u = vec(m);
    code = SVDGetSingularTriplet(svd, i, &svalue, u, v); CHKERR(code);
    MatSetValues(V, n, &range(n)[0], 1, &i, &stl(v)[0], INSERT_VALUES);
    MatSetValues(U, m, &range(m)[0], 1, &i, &stl(u)[0], INSERT_VALUES);
    MatSetValue(S, i, i, 1.0/svalue, INSERT_VALUES);
  }
  MatAssemblyBegin(V, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(V, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(U, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(U, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY);
  SVDDestroy(&svd); 

  /* V S^{-1} U^T */
  Mat Ut = mat(l, n);
  code = MatTranspose(U, MAT_INITIAL_MATRIX, &Ut); CHKERR(code);

  Mat B = mat(l, n, "dense");
  MatMatMult(S, Ut, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &B);

  Mat C = mat(n, n, "dense");
  MatMatMult(V, B, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C);

  return C;
}

Vec projectOn(Vec u, vector<Vec> v){
  for (int i = 0; i < v.size(); i++){

  }
}
