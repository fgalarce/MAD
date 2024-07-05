#include <iostream>
#include <math.h>
#include <vector>
#include <slepc.h>
#include <STLvectorUtils.hpp>
#include <slepcWrapper.hpp>
#include <tools.hpp>
#include <assert.h>
#include <parameters.hpp>
#include <io.hpp>

using namespace std;

int main(int argc, char *argv[], char *envp[]){

  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL);

  /* Parse input data */
  assert(argc == 4);
  Parameters par(argv[3]);
  par.print();
  int coord_i = stoi(argv[1]); /* time coordinate */
  int coord_j = stoi(argv[2]); /* heart rate coordinate */

  /* Initialize MDMA */
  IO io;
  io.initialize(par);
  int nbDofs = io.nbVertices()*3;

  /* Get heart Rate of the patients */
  vector<double> heartRateData(par.nbSimulations());
  for (size_t i = 0; i < par.nbSimulations(); i++){
    heartRateData[i] = stod(getParameter("heartRate", par.maniFolder() + "/sim" + wildcard(i) + "/data"));
    cout << "HR sim " << i << ": " << heartRateData[i] << endl;
  }
  double HR_max = max(heartRateData);
  double HR_min = min(heartRateData);

  /* Gather folders with all simulations to be loaded */
  double HR_window_left = HR_min + coord_j*(HR_max - HR_min)/par.nbHeartRateWindows();
  double HR_window_right = HR_min + (coord_j+1)*(HR_max - HR_min)/par.nbHeartRateWindows();

  /* Search for simulations in heart rate range*/
  vector<int> simsToLoad;
  for (int i = 0; i < par.nbSimulations(); i++){
    if (heartRateData[i] <= HR_window_right && heartRateData[i] >= HR_window_left ){
      simsToLoad.push_back(i);
    }
  }


  /* Compute number of snapshots. This could look cumbersome but saves RAM memory because snapshot matrix can be pre-allocated */
  vector<vector<size_t>> timeIndexes(simsToLoad.size());
  size_t nbSnapshots = 0;
  for (int i = 0; i < simsToLoad.size(); i++){
    
    double t_min = 1.0/double(par.nbTimeWindows())*(coord_i)*60/heartRateData[simsToLoad[i]];
    double t_max = 1.0/double(par.nbTimeWindows())*(coord_i+1)*60/heartRateData[simsToLoad[i]];
    double time = t_min;
    while (time < t_max && time < ((0.99)*60/heartRateData[simsToLoad[i]]) ){
      timeIndexes[i].push_back(int(floor(time/par.timeStep())));
      time = time + par.timeStep();
    }
    nbSnapshots = nbSnapshots + timeIndexes[i].size();
  }

  /* Snapshots */
  Mat A = mat(nbDofs, nbSnapshots, "dense"); 
  Vec V = vec(nbDofs);

  PetscInt iterSnap = 0;
  PetscInt iterSim = 0;

  PetscErrorCode code;
  for (auto it = simsToLoad.cbegin(); it!= simsToLoad.cend(); it++){
    for (auto itTime = timeIndexes[iterSim].cbegin(); itTime != timeIndexes[iterSim].cend(); itTime++){
      loadVec(V, par.maniFolder() + "sim" + wildcard(*it) + "/velocity." + wildcard(*itTime) + par.inputFormat());
      double snapshot[nbDofs];
      code = VecGetValues(V, nbDofs, &range(nbDofs)[0], snapshot); CHKERRQ(code); 
      MatSetValues(A, nbDofs, &range(nbDofs)[0], 1, &iterSnap, snapshot, INSERT_VALUES);
      iterSnap = iterSnap + 1;
    }
    iterSim = iterSim + 1;
  }
  VecDestroy(&V);

  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  cout << "\n----------------------------------" << endl;
  cout << "Creating POD basis for window (" << coord_i << ", " << coord_j << ")" << endl;
  cout << "Time (normalized in [0,1]) : [" << 1.0/double(par.nbTimeWindows())*coord_i << ", " << 1.0/double(par.nbTimeWindows())*(coord_i+1) << "] seg" << endl;
  cout << "Heart rate : [" << HR_window_left << ", " << HR_window_right << "]" << endl;
  cout << "HRmax: " << HR_max << ". HRmin: " << HR_min << endl;
  cout << "# snapshots : " << nbSnapshots << " (" << simsToLoad.size() << " simulation(s))" << endl;

  PetscScalar normSnapshots;
  MatNorm(A, NORM_FROBENIUS, &normSnapshots);
  cout << "Norm snapshots matrix: " << normSnapshots << endl;

  /* Compute co-variance matrix */
  Mat A_tr = mat(nbSnapshots, nbDofs);
  code = MatTranspose(A, MAT_INITIAL_MATRIX, &A_tr); CHKERRQ(code);
  MatAssemblyBegin(A_tr, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A_tr, MAT_FINAL_ASSEMBLY);

  Mat covMatrix = mat(nbSnapshots, nbSnapshots);
  Mat covMatrix_aux;

  Mat massMatrix, stiffnessMatrix;
  if (par.innerProduct() == "l2"){
    code = MatMatMult(A_tr, A, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &covMatrix); 
  } else if (par.innerProduct() == "L2"){

    /* Import mass matrix */
    mat(massMatrix, nbDofs, nbDofs);
    loadMat(massMatrix, par.mass_matrix_path());

    mat(covMatrix_aux, nbSnapshots, nbDofs);
    code = MatMatMult(A_tr, massMatrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &covMatrix_aux); 
    code = MatMatMult(covMatrix_aux, A, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &covMatrix); 

  } else if (par.innerProduct() == "H1"){

    /* Import mass matrix and stifness matrix */
    mat(massMatrix, nbDofs, nbDofs);
    mat(stiffnessMatrix, nbDofs, nbDofs);
    loadMat(massMatrix, par.mass_matrix_path());
    loadMat(stiffnessMatrix, par.stiffness_matrix_path());

    code = MatAXPY(massMatrix, 1.0, stiffnessMatrix, SAME_NONZERO_PATTERN); /* M + K */

    cout << "Computing A^T * (K + M) * A" << endl;
    code = MatMatMult(A_tr, massMatrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &covMatrix_aux); 
    code = MatMatMult(covMatrix_aux, A, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &covMatrix); 

    PetscReal normMatrix;
    MatNorm(massMatrix, NORM_FROBENIUS, &normMatrix);
    cout << "Norm mass: " << normMatrix << endl;
    MatNorm(stiffnessMatrix, NORM_FROBENIUS, &normMatrix);
    cout << "Norm stiffness: " << normMatrix << endl;
  }
  
  MatAssemblyBegin(covMatrix, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(covMatrix, MAT_FINAL_ASSEMBLY);

  MatNorm(covMatrix, NORM_FROBENIUS, &normSnapshots);
  cout << "Norm covariance matrix: " << normSnapshots << endl;

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
  cout << "nconv : " << nconv << endl; 

  for (int i=0; i<par.nbModes(); i++) {
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

    if (par.outputFormat() == ".ens" || par.outputFormat() == ".vct"){
      io.writeState(mode, "mode", i);
    } else if (par.outputFormat() == ".bin"){ 
      saveVec(mode, fileNameOut + ".bin");
    }
    code = VecDestroy(&mode);
    code = VecDestroy(&Vi); CHKERRQ(code);
    code = VecDestroy(&Vr); CHKERRQ(code);
  }

  string fileNameData = par.dirResults() + "/snap_data.txt";
  ofstream outFile(fileNameData);
  cout << "Writing: " << fileNameData << endl;
  outFile << "--- LOG window id: " << to_string(coord_i) + "_" + to_string(coord_j) << " ---\n";
  outFile << "Number of Snapshots: " << nbSnapshots << endl;
  outFile << "Simulation Id | HR | time : " << endl;
  for (int simId = 0; simId < simsToLoad.size(); simId ++){
    outFile << simsToLoad[simId] << " " << heartRateData[simId];
    for (int timeId = 0; timeId < timeIndexes[simId].size() ; timeId++){
      outFile << " " << timeIndexes[simId][timeId];
    }
    outFile << endl;
  }
  outFile.close();

  exportData(par.dirResults() + "./eigenValues.txt", eValues);

  if (par.innerProduct() == "L2" || par.innerProduct() == "H1"){
    code = MatDestroy(&massMatrix); CHKERRQ(code); 
    code = MatDestroy(&covMatrix_aux); CHKERRQ(code); 
  }

  if (par.innerProduct() == "H1")
    code = MatDestroy(&stiffnessMatrix); CHKERRQ(code); 

  MatDestroy(&A);
  MatDestroy(&A_tr);  
  MatDestroy(&covMatrix);  
  EPSDestroy(&eps);
  SlepcFinalize();
}
