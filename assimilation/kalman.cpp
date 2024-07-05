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

#include<kalman.hpp>

void Kalman::initialize(Parameters parameters, const Geometry & geometry, Model model, const InnerProduct & innerProduct, const Measures usMeasures){
  nbDofs = geometry.nbVertices * 3;
  geo = geometry;
  Vn = model;
  par = parameters;
  m_phi = model.basis(0.0);
  m_ip = innerProduct;
  m_nbModes = par.nbModes();
  m_nbMeasures = usMeasures.nbMeasures();
  
}

void Kalman::initialize(Parameters parameters, const Geometry & geometry, const InnerProduct & innerProduct, const Measures usMeasures){
  nbDofs = geometry.nbVertices * 3;
  m_nbMeasures = usMeasures.rieszRepresenters().size();
}

Vec Kalman::mean(const vector<Vec> & sample){
  int vecSize;
  VecGetSize(sample[0], &vecSize);
  Vec average = zeros(vecSize);
  for (int i=0;i<sample.size();i++){
    VecAXPY(average, 1.0/sample.size(), sample[i]); 
  }
  return average;
}

Mat Kalman::covariance(const vector<Vec> & sample){
  int n;
  VecGetSize(sample[0], &n);
  Mat C = mat(n, n, "dense");
  Vec average = mean(sample);
  for (int i=0;i<sample.size();i++){
    Mat C_i = outer(sample[i]);
    MatAXPY(C, 1.0, C_i, DIFFERENT_NONZERO_PATTERN);
  }
  MatScale(C, 1.0/(sample.size() - 1));
  return C;
}

vector<Vec> Kalman::samples(int iteration){
  vector<Vec> a(par.nbSamples());
  vector<Vec> samples(par.nbSamples());
  for (int i=0; i<par.nbSamples(); i++){
    samples[i] = vec(nbDofs, par.maniFolder() + "/sim" + wildcard(i) + "/stokes_vector." + wildcard(iteration) + ".vct");
    /* Project snapshot on V_n */
    vec(a[i], par.nbModes());
    for (int j=0; j<par.nbModes(); j++){
      double proj_Vn_u = m_ip(m_phi[j], samples[i]);
      VecSetValues(a[i], 1, &j, &proj_Vn_u, INSERT_VALUES);
    }
    VecAssemblyBegin(a[i]);
    VecAssemblyEnd(a[i]);
  }
  return a;
}

/* Compute A C A^T */
Mat Kalman::propagateCovariance(Mat C, Mat A){

  int m,n;
  MatGetSize(A, &m, &n);
  assert(m == n and "Covariance not square");
  PetscErrorCode code;

  Mat D = mat(n, n);
  MatMatMult(A, C, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &D);

  Mat Atr = mat(n, m);
  code = MatTranspose(A, MAT_INITIAL_MATRIX, &Atr); CHKERR(code);

  Mat E = mat(n, n);
  MatMatMult(D, Atr, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &E);

  return E;
}

Mat Kalman::assembleFilter(Mat C10, Mat G, Mat Gt, Mat T){
  Mat M = mat(m_nbMeasures, m_nbModes);
  MatMatMult(G, C10, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &M);

  Mat N = mat(m_nbMeasures, m_nbMeasures);
  MatMatMult(M, Gt, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &N);
  MatAXPY(N, 1.0, T, DIFFERENT_NONZERO_PATTERN);

  Mat O = mat(m_nbModes, m_nbModes);
  MatCopy(C10, O, DIFFERENT_NONZERO_PATTERN);

  Mat P = mat(m_nbModes, m_nbMeasures);
  MatMatMult(O, Gt, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &P);

  Mat K = mat(m_nbModes, m_nbMeasures);
  MatMatMult(P, inverse(N), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &K);

  return K;
}
