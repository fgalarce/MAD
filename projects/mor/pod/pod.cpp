/*=============================================================================
  This file is part of the code MAD
  Copyright (C) 2017-2024,
    
     Felipe Galarce

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

#include <mad.hpp>

int main(int argc, char *argv[], char *envp[]){

  /* Parse input data */
  Parameters par = MADinitialize(argc, argv); 

  PetscErrorCode code;

  /* Initialize MDMA */
  IO io;
  io.initialize(par);

  Geometry geo;
  geo.initialize(par, io);

  Boundary bd;
  bd.initialize(par, geo);

  InnerProduct ip; 
  ip.initialize(par, geo, bd);

  /* Computing number of snaps */
  int jump = 1;
  int nbSnapshots = 0;
  bool static_problem = false;
  if (static_problem){
    nbSnapshots = par.nbSimulations();
  } else {
    for (int simId = 0; simId < par.nbSimulations(); simId++){
      if (par.snaps_per_sim() == -1){
        nbSnapshots = nbSnapshots + floor(stoi(getParameter("nbIterations", par.maniFolder() + "sim" + wildcard(simId) + "/parameters"))/jump) - par.start();
      } else {
        nbSnapshots = nbSnapshots + par.snaps_per_sim();
      }
    }
  }

  /* Snapshots */
  int nbDofs = 0;
  for (int i = 0; i < par.nbVariables(); i++){
    nbDofs = nbDofs + io.nbVertices()*par.nbDofsPerNode()[i];
  }
  vector<Vec> V(nbSnapshots);
  int m,n;
  Mat A = mat(nbDofs, nbSnapshots, "dense");
  code = MatGetOwnershipRange(A, &m, &n); CHKERR(code);
  int nbDofsLocal = n - m;

  MADprint("POD: Assembling matrix with # snapshots : ", nbSnapshots); 
  int iterSnap = 0;
  for (int simId = 0; simId < par.nbSimulations(); simId++){
    int nbSnaps;
    if (static_problem){
      nbSnaps = 1;
    } else {
      if (par.snaps_per_sim() == -1){
        nbSnaps = floor(stoi(getParameter("nbIterations", par.maniFolder() + "sim" + wildcard(simId) + "/parameters")));
      } else {
        nbSnaps = par.snaps_per_sim() + par.start();
      }
    }
    for (int i = par.start(); i < nbSnaps; i = i + jump){
      vec(V[iterSnap], nbDofs);
      vector<string> to_load(par.nbVariables());
      for (int var = 0; var < par.nbVariables(); var++){
        string end_str;
        if (par.nbDofsPerNode()[var] == 1){
          to_load[var] = par.maniFolder() + "sim" + wildcard(simId) + "/" + par.variableName()[var] + "." + wildcard(i) + ".scl"; 
        } else {
          to_load[var] = par.maniFolder() + "sim" + wildcard(simId) + "/" + par.variableName()[var] + "." + wildcard(i) + ".vct"; 
        }
      }
      io.loadVector(V[iterSnap], to_load);
      code = VecAssemblyBegin(V[iterSnap]); CHKERRQ(code);
      code = VecAssemblyEnd(V[iterSnap]); CHKERRQ(code);
      double snapshot[nbDofsLocal];
      code = VecGetValues(V[iterSnap], nbDofsLocal, &range(m, n)[0], snapshot); CHKERRQ(code); 
      code = MatSetValues(A, nbDofsLocal, &range(m, n)[0], 1, &iterSnap, snapshot, INSERT_VALUES); CHKERRQ(code);

      iterSnap = iterSnap + 1;
    }
  }
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  MADprint("POD: Total number of snapshots: ", iterSnap);

  double normSnapshots;
  code = MatNorm(A, NORM_FROBENIUS, &normSnapshots); CHKERRQ(code);
  MADprint("POD: Norm snapshots matrix: ", normSnapshots);
  MADprint("POD: Creating POD basis.");
  MADprint("POD: assembling A^T M A");
  Mat covMatrix = mat(nbSnapshots, nbSnapshots, "dense");
  int cov_m, cov_n;
  code = MatGetOwnershipRange(covMatrix, &cov_m, &cov_n); CHKERRQ(code);
  for (int i = 0; i < nbSnapshots; i++){
    for (int j = i; j < nbSnapshots; j++){
      double iprod = ip(V[i], V[j]);
      if (i >= cov_m && i < cov_n){
        if(i == j){
          iprod = iprod / 2.0; // this is necessary because the loop pass by here two times
        }
        code = MatSetValue(covMatrix, i, j, iprod, INSERT_VALUES); CHKERRQ(code);
        cout << "     row-col: " << scientific << i + 1 << " - " << j + 1 << " / " << nbSnapshots << ". ip = " << iprod << "   \r";
      }
    }
  }
  MADprint("\n");
  code = MatAssemblyBegin(covMatrix, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyEnd(covMatrix, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

  /* Only upper triangular matrix was assembled, so C = 0.5 ( C + C^T) */
  Mat covMatrix_tr = mat(nbSnapshots, nbSnapshots, "dense");
  code = MatTranspose(covMatrix, MAT_INITIAL_MATRIX, &covMatrix_tr); CHKERRQ(code);
  code = MatAssemblyBegin(covMatrix_tr, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyEnd(covMatrix_tr, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAXPY(covMatrix, 1.0, covMatrix_tr, DIFFERENT_NONZERO_PATTERN); CHKERRQ(code);

  double normCov;
  code = MatNorm(covMatrix, NORM_FROBENIUS, &normCov); CHKERRQ(code);
  MADprint("POD: norm covariance matrix: ", normCov);

  /* Compute eigen-vectors */
  vector<double> eValues, eValues_i;
  eValues.resize(par.nbModes());
  eValues_i.resize(par.nbModes());

  EPS eps;
  PetscInt nconv;
  EPSCreate(PETSC_COMM_WORLD, &eps);
  EPSSetOperators(eps, covMatrix, NULL);
  EPSSetProblemType(eps, EPS_GHEP);
  EPSSetType(eps, EPSKRYLOVSCHUR);
  EPSSetDimensions(eps, par.nbModes(), PETSC_DECIDE, PETSC_DECIDE);
  EPSSetFromOptions(eps);
  EPSSolve(eps);
  EPSGetConverged(eps, &nconv);
  MADprint("nconv : ", nconv); 

  vector<string> modeNames(par.nbVariables());
  for (int i = 0; i < par.nbVariables(); i++){
    modeNames[i] = par.variableName()[i] + "_mode";
  }

  for (int i=0; i<par.nbModes(); i++){
    Vec Vr = vec(nbSnapshots);
    Vec Vi = vec(nbSnapshots);
    code = EPSGetEigenvalue(eps, i, &eValues[i], &eValues_i[i]); CHKERRQ(code);
    code = EPSGetEigenvector(eps, i, Vr, Vi); CHKERRQ(code);
  
    string fileNameOut = par.dirResults() + "/mode." + wildcard(i);

    /* Compute mode */
    Vec mode = vec(nbDofs); 
    code = MatMult(A, Vr, mode); CHKERRQ(code);
    PetscScalar singValue;
    singValue = sqrt(eValues[i]);
    code = VecScale(mode, 1.0/singValue); CHKERRQ(code);

    io.writeState(mode, i, modeNames);

    code = VecDestroy(&mode);
    code = VecDestroy(&Vi); CHKERRQ(code);
    code = VecDestroy(&Vr); CHKERRQ(code);
  }

  exportData(par.dirResults() + "./eigenValues.txt", eValues);

  MatDestroy(&covMatrix);  
  EPSDestroy(&eps);
  MADfinalize(par);
}
