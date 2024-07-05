#include <ultra-4d-flow.hpp>

/* 
 * Receives a non divergence free field and project it by using a potential \tsi such that
 * u_out = \grad \tsi + u_in
 * It requires a scalar stiffnes matrix to solve an homogeneous Laplacian
 * */

int main(int argc, char *argv[]){

  assert(argc == 2);

  string data_file = argv[1];
  
  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  Parameters par(data_file);
  par.print();

  IO io;
  io.initialize(par);

  InnerProduct ip;

  int nbVertices = io.nbVertices();
  int nbDofs = 3 * io.nbVertices(); 

  Vec u_in = vec(nbDofs);
  Vec u_out = vec(nbDofs);
  loadVec(u_in, par.dirSyntheticField() + "/test1_caro42kv/up_vector.vct");
//  loadVec(u_out, par.dirSyntheticField() + "/test2_caro42kv/up_vector.vct");
  loadVec(u_out, par.dirSyntheticField() + "/test3_caro42kv/up_vector.vct");

  Vec v = vec(nbDofs);
  code = VecCopy(u_in, v); CHKERRQ(code);
  code = VecAXPY(v, -1.0, u_out); CHKERRQ(code);


  /* Import Riesz representers */
  int m = 216;
  Mat W = mat(nbDofs, m, "dense");
  for (int i = 0; i < m; i++){
    Vec rr = vec(nbDofs);
    loadVec(rr, par.rieszRepresentersFolder() + "/rieszRepresenter." + wildcard(i) + ".bin");

    vector<double> rr_vec(nbDofs);
    VecGetValues(rr, nbDofs, &range(nbDofs)[0], &rr_vec[0]);
    MatSetValues(W, nbDofs, &range(nbDofs)[0], 1, &i, &rr_vec[0], INSERT_VALUES);

    VecDestroy(&rr);
  }
  code = MatAssemblyBegin(W, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyEnd(W, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  Mat W_tr = mat(m, nbDofs);
  code = MatTranspose(W, MAT_INITIAL_MATRIX, &W_tr); CHKERRQ(code);

  /* W^T (M + K) W*/
  Mat F = mat(nbDofs, nbDofs);
  loadMat(F, par.mass_matrix_path());

  Mat K = mat(nbDofs, nbDofs);
  loadMat(K, par.stiffness_matrix_path());

  code = MatAXPY(F, 1.0, K, SAME_NONZERO_PATTERN); CHKERRQ(code);

  Mat Aux = mat(nbDofs, m);
  code = MatMatMult(F, W, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Aux); CHKERRQ(code);

  Mat A = mat(m, m);
  code = MatMatMult(W_tr, Aux, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &A); CHKERRQ(code);

  KSP ksp = configureKSP(A, "preonly", "lu");
  Vec c = vec(m);
  Vec aux = vec(nbDofs);
  code = MatMult(F, v, aux); CHKERRQ(code);
  Vec rhs = vec(m);
  code = MatMult(W_tr, aux, rhs); CHKERRQ(code);

  cout << "Norm_fro of LHS: " << norm(A) << endl;
  cout << "Norm_2 of RHS: " << norm(rhs) << endl;
  
  code = KSPSolve(ksp, rhs, c); CHKERRQ(code); 

  Vec error = vec(nbDofs);
  code = MatMult(W, c, error); CHKERRQ(code);

  cout << norm(c) << endl;

  VecAXPY(error, -1.0, v);

  io.writeState(v, "rr_vitesse");

  cout << "Norm residual: " << sqrt(ip(error, F, error)) << endl;
  cout << "Norm RR(u,p): " << sqrt(ip(v,F,v)) << endl;
  SlepcFinalize();


}
