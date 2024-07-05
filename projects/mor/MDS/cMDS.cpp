#include<ultra-4d-flow.hpp>

int main(int argc, char *argv[]){

  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  /* Parse data file */
  assert(argc == 3);
  Parameters par(argv[1]);

  /* Initialize MDMA objects */
  par.print();

  int nbMeshes = par.nbMeshes();
  cout << "Assembling dissimilarity matrix" << endl;
  Mat D = mat(nbMeshes, nbMeshes, "dense");
  for (int i = 0; i < nbMeshes; i++){
    for (int j = 0; j < nbMeshes; j++){
      if (i != j){
        vector<double> data0 = importdata1D(par.dirResults() + "/D" + to_string(i) + "_" + to_string(j) + ".txt");
        vector<double> data1 = importdata1D(par.dirResults() + "/D" + to_string(j) + "_" + to_string(i) + ".txt");
        code = MatSetValue(D, i, j, (0.5*(data0[0] + data1[0])*0.5*(data0[0] + data1[0])), INSERT_VALUES); CHKERRQ(code);
      } else {
        code = MatSetValue(D, i, j, 0.0, INSERT_VALUES); CHKERRQ(code);
      }
    }
  }
  code = MatAssemblyBegin(D, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyEnd(D, MAT_FINAL_ASSEMBLY); CHKERRQ(code);


  cout << "Assembling X'X" << endl;
  /* Assemble H = I - 1/nbMeshes e e^T*/
  Mat H = eye(nbMeshes, "dense");
  Vec uno = ones(nbMeshes);
  code = MatAXPY(H, -1.0/nbMeshes, outer(uno), DIFFERENT_NONZERO_PATTERN); CHKERRQ(code);

  VecView(uno, PETSC_VIEWER_STDOUT_SELF);
  MatView(outer(uno), PETSC_VIEWER_STDOUT_SELF);
  MatView(H, PETSC_VIEWER_STDOUT_SELF);

  /* Assemble X'X = HDH'*/
  Mat XX = mat(nbMeshes, nbMeshes, "dense");
  Mat XX0 = mat(nbMeshes, nbMeshes, "dense");
  MatMatMult(H, D, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &XX0);

  Mat Htr = mat(nbMeshes, nbMeshes, "dense");
  code = MatTranspose(H, MAT_INITIAL_MATRIX, &Htr); CHKERR(code);
  code = MatMatMult(XX0, Htr, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &XX); CHKERRQ(code);
  code = MatScale(XX, -1.0/2.0); CHKERRQ(code);

  /* Compute cMDS */
  vector<double> eValues, eValues_i, svalues;
  eValues.resize(nbMeshes);
  eValues_i.resize(nbMeshes);
  svalues.resize(nbMeshes);
  vector<Vec> modes(nbMeshes);

  EPS eps;
  PetscInt nconv;
  EPSCreate(PETSC_COMM_WORLD, &eps);
  EPSSetOperators(eps, XX, NULL);
  EPSSetProblemType(eps, EPS_HEP);
  EPSSetType(eps, EPSKRYLOVSCHUR);
  EPSSetDimensions(eps, nbMeshes, PETSC_DECIDE, PETSC_DECIDE);
  EPSSetFromOptions(eps);
  EPSSolve(eps);
  EPSGetConverged(eps, &nconv);
  cout << "nconv : " << nconv << endl; 

  for (int i=0; i<nbMeshes; i++){
    Vec Vi = vec(nbMeshes);
    vec(modes[i], nbMeshes);
    code = EPSGetEigenvalue(eps, i, &eValues[i], &eValues_i[i]); CHKERRQ(code);
    code = EPSGetEigenvector(eps, i, modes[i], Vi); CHKERRQ(code);
  
    /* Compute mode */
    if (eValues[i] > 0){
      svalues[i] = sqrt(eValues[i]);
    } else {
      svalues[i] = 0.0;
    }
    cout << eValues[i] << endl;
//    code = VecScale(modes[i], singValue); CHKERRQ(code);

    code = VecDestroy(&Vi); CHKERRQ(code);
  }

  cout << "Computing low-dimensional representation." << endl;
  int dimension = par.nbModes();
  vector<Vec> Y(nbMeshes);
  vector<vector<double>> mode_low(nbMeshes);
  for (int i = 0; i < nbMeshes; i++){
    vec(Y[i], dimension);
    mode_low[i].resize(dimension);
    VecGetValues(modes[i], dimension, &range(dimension)[0], &mode_low[i][0]);
    if (eValues[i] > 0){
      mode_low[i] = svalues[i]*mode_low[i];
    } else {
      cout << "Warning: non-positive eigen-value: ";
      cout << scientific << eValues[i] << endl;
    }
    VecSetValues(Y[i], dimension, &range(dimension)[0], &mode_low[i][0] , INSERT_VALUES);
  }
  exportData("./Y_work.txt", mode_low, " ",3);

  cout << "Computing Euclidean distances." << endl;
  vector<vector<double>> d_ij(nbMeshes);
  for (int i = 0; i < nbMeshes; i++){
    d_ij[i].resize(nbMeshes);
    for (int j = 0; j < nbMeshes; j++){
      Vec dist = vec(dimension);
      VecCopy(Y[i], dist);
      VecAXPY(dist, -1.0, Y[j]);
      d_ij[i][j] = norm(dist);
    }
  }

  exportData(par.dirResults() + "/Dhat_work" + argv[2] + ".txt", d_ij, " ", 3);
  
  EPSDestroy(&eps);
  par.finalize();
  SlepcFinalize();
}
