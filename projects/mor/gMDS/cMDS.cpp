#include<ultra-4d-flow.hpp>

int main(int argc, char *argv[]){

  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  /* Parse data file */
  assert(argc == 2);
  Parameters par(argv[1]);

  /* Initialize MDMA objects */
  par.print();

  /* Dissimilarity matrix */
  Mat D = mat(par.nbMeshes(), par.nbMeshes(), "dense");

  IO io;
  io.initialize(par);

  /* Load surface registration mappings */
  for (int i = 0; i < par.nbMeshes(); i++){

    IO io_i;
    io_i.initialize(par.geometryData() + "/venturi" + wildcard(i) + ".mesh");

    Geometry geo;
    geo.initialize(par, io_i);

    int nbDofs = 3*geo.nbVertices;

    /* Compute work done by LDDMM field */
    CFD cfd;
    cfd.initialize(par, geo); 

    /* T^T Mb T */
    cfd.assembleMassBoundary();

    for (int j = i; j < par.nbMeshes(); j++){
  
      if (i != j){
        Vec lddmm_kl = vec(nbDofs);
        loadVec(lddmm_kl, par.surfaceMapping() + "/" + wildcard(i) + "/lddmm_" + wildcard(j) + ".vct");

        Vec lddmm_lk = vec(nbDofs);
        loadVec(lddmm_lk, par.surfaceMappingInverse() + "/mesh" + wildcard(j) + "/" + wildcard(i) + "/I_Tij.00000.vct");

        Vec mtv = vec(nbDofs);
        code = MatMult(cfd.massBD(), lddmm_kl, mtv); CHKERRQ(code);
        double dij1;
        code = VecDot(mtv, lddmm_kl, &dij1); CHKERRQ(code);

        mtv = zeros(nbDofs);
        code = MatMult(cfd.massBD(), lddmm_lk, mtv); CHKERRQ(code);
        double dij2;
        code = VecDot(mtv, lddmm_lk, &dij2); CHKERRQ(code);

        double dij = (sqrt(dij1) + sqrt(dij2))/2;
        dij = dij*dij;
      
        code = MatSetValue(D, i, j, dij, INSERT_VALUES); CHKERRQ(code);
        code = MatSetValue(D, j, i, dij, INSERT_VALUES); CHKERRQ(code);
      } else {
        code = MatSetValue(D, i, j, 0.0, INSERT_VALUES); CHKERRQ(code);
      }
    }
  }

  code = MatAssemblyBegin(D, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyEnd(D, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  saveMat(D, "D.m");


  cout << "Assembling X'X" << endl;
  /* Assemble H = I - 1/nbMeshes e e^T*/
  Mat H = eye(par.nbMeshes(), "dense");
  Vec uno = ones(par.nbMeshes());
  code = MatAXPY(H, -1.0/par.nbMeshes(), outer(uno), DIFFERENT_NONZERO_PATTERN); CHKERRQ(code);

  /* Assemble X'X = HDH */
  Mat XX = mat(par.nbMeshes(), par.nbMeshes(), "dense");
  Mat XX0 = mat(par.nbMeshes(), par.nbMeshes(), "dense");
  code = MatMatMult(H, D, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &XX0);
  code = MatMatMult(XX0, H, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &XX); CHKERRQ(code);
  code = MatScale(XX, -1.0/2.0); CHKERRQ(code);

  cout << "Computing low-dimensional representation." << endl;
  /* Compute cMDS */
  vector<double> eValues(par.nbMeshes());
  vector<double> svalues;
  vector<vector<double>> modes; 

  EPS eps;
  int nconv;
  EPSCreate(PETSC_COMM_WORLD, &eps);
  EPSSetOperators(eps, XX, NULL);
  EPSSetProblemType(eps, EPS_HEP);
  EPSSetType(eps, EPSKRYLOVSCHUR);
  EPSSetDimensions(eps, par.nbMeshes(), PETSC_DECIDE, PETSC_DECIDE);
  EPSSetFromOptions(eps);
  EPSSolve(eps);
  EPSGetConverged(eps, &nconv);
  cout << "EPS: nconv : " << nconv << endl; 

  for (int i=0; i<par.nbMeshes(); i++){
    code = EPSGetEigenvalue(eps, i, &eValues[i], NULL); CHKERRQ(code);
  
    /* Compute mode */
    if (eValues[i] > 0){
      Vec Vr = vec(par.nbMeshes());
      code = EPSGetEigenvector(eps, i, Vr, NULL); CHKERRQ(code);
      modes.push_back(stl(Vr));
      svalues.push_back(sqrt(eValues[i]));
      cout << "gMDS: mode" << i << ". EV = " << eValues[i] << endl;
      cout << modes.back() << endl;
    }
  }

  for (int i=0; i<svalues.size(); i++){
    modes[i] = svalues[i] * modes[i];
  }

  int dimension = svalues.size();

  cout << "gMDS: Maximal Euclidean dimension: " << dimension << endl;

  vector<vector<double>> Dij_hat(par.nbMeshes());


  for (int i = 0; i < par.nbMeshes(); i++){
    Dij_hat[i].resize(par.nbMeshes());
    for (int j = 0; j < par.nbMeshes(); j++){
      Dij_hat[i][j] = norm(getCol(modes, i) - getCol(modes, j));
      Dij_hat[i][j] = Dij_hat[i][j] * Dij_hat[i][j];
    }
  }

  double Dfro, Dfro_error;
  MatNorm(D, NORM_FROBENIUS, &Dfro);
  Mat error = mat(par.nbMeshes(), par.nbMeshes(), "dense");
  MatCopy(D, error, DIFFERENT_NONZERO_PATTERN);
  MatAXPY(error, -1.0, petsc(Dij_hat), DIFFERENT_NONZERO_PATTERN);
  MatNorm(error, NORM_FROBENIUS, &Dfro_error); 

  cout << "cMDS: Fro error = " << Dfro_error / Dfro;

  exportData("D_hat.txt", Dij_hat);
 
  EPSDestroy(&eps);

  par.finalize();
  SlepcFinalize();
}
