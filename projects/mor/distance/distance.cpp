#include<ultra-4d-flow.hpp>

int main(int argc, char *argv[]){

  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  /* Parse data file */
  assert(argc == 2);
  Parameters par(argv[1]);

  /* Initialize MDMA objects */
  par.print();

  IO io, io_target;
  io.initialize(par);

  Geometry geo;
  geo.initialize(par, io);

  InnerProduct ip;
  ip.initialize(par, geo);

  LinearAlgebra la;
  la.initialize(par, ip);

  /* Load POD basis */
  int n = par.nbModes();
  int nbDofs = geo.nbVertices * 3;
  vector<Vec> phi0(n);
  vector<Vec> phi1(n);
  for(int i=0;i<n;i++){
    vec(phi0[i], nbDofs);
    vec(phi1[i], nbDofs);
    loadVec(phi0[i], par.templateModel() + "/mode." + wildcard(i)+".vct");
    loadVec(phi1[i], par.model() + "/pod." + wildcard(i) + ".vct");
  }
  la.orthonormalize(phi0);
  la.orthonormalize(phi1);

  if (!la.checkOrthonormality(phi0)){
    cout << "MAIN: basis from " << par.templateModel() << " are not orthonormal." << endl;
    exit(1);
  } 
  if (!la.checkOrthonormality(phi1)){
    cout << "MAIN: basis from " << par.model() << " are not orthonormal." << endl;
    exit(1);
  }

  Mat Phi0 = buildProjector(phi0);
  Mat Phi1 = buildProjector(phi1);
  Mat Phi0_tr = mat(n, nbDofs, "dense");
  Mat Phi1_tr = mat(n, nbDofs, "dense");

  cout << "MAIN: Norm P_0: " << norm(Phi0) << endl;
  cout << "MAIN: Norm P_1: " << norm(Phi1) << endl;

  cout << "MAIN: Solving GEVP to compute distance" << endl;
  code = MatTranspose(Phi0, MAT_INITIAL_MATRIX, &Phi0_tr); CHKERRQ(code);
  code = MatTranspose(Phi1, MAT_INITIAL_MATRIX, &Phi1_tr); CHKERRQ(code);

  Mat Phi0_tr_M = mat(n, nbDofs, "dense");
  code = MatMatMult(Phi0_tr, ip.matrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Phi0_tr_M); CHKERRQ(code);

  Mat P_Phi0_M = mat(nbDofs, nbDofs, "dense");
  code = MatMatMult(Phi0, Phi0_tr_M, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &P_Phi0_M); CHKERRQ(code);
  
  Mat Phi1_tr_M = mat(n, nbDofs, "dense");
  code = MatMatMult(Phi1_tr, ip.matrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Phi1_tr_M); CHKERRQ(code);

  Mat P_Phi1_M = mat(nbDofs, nbDofs, "dense");
  code = MatMatMult(Phi1, Phi1_tr_M, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &P_Phi1_M); CHKERRQ(code);

  code = MatAXPY(P_Phi1_M, -1.0, P_Phi0_M, DIFFERENT_NONZERO_PATTERN); CHKERRQ(code);

  Mat A = mat(nbDofs, nbDofs, "dense");
  code = MatMatMult(ip.matrix, P_Phi1_M, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &A);
  
  Mat B = mat(nbDofs, nbDofs, "dense");
  code = MatMatMult(A, A, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &B);

  cout << "MAIN: Norm M: " << norm(ip.matrix) << endl;
  vector<vector<double>> e_values = la.eigenValueProblem(B, ip.matrix, "GHEP"); 

  cout << "MAIN: d(E,F) = " << max(e_values[0]) << endl; 

  ofstream file(par.dirResults() + "/distance_ij.txt");
  file << max(e_values[0]) << endl;
  
  
  MatDestroy(&Phi1);
  MatDestroy(&Phi0);
  MatDestroy(&Phi0_tr);
  MatDestroy(&Phi1_tr);
  MatDestroy(&A);
  SlepcFinalize();
}
