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

#include <petscWrapper.hpp>

using namespace std;

void matSet(Mat A, int i, int j, double Aij){
  int m,n;
  PetscErrorCode code;
  MatGetOwnershipRange(A, &m, &n); 
  if (i >= m && i < n){
    code = MatSetValue(A, i, j, Aij, ADD_VALUES); CHKERR(code);
  }
}

void matSetInsert(Mat A, int i, int j, double Aij){
  int m,n;
  MatGetOwnershipRange(A, &m, &n); 
  if (i >= m && i < n){
    MatSetValue(A, i, j, Aij, INSERT_VALUES);
  }
}

void vecSet(Vec u, int i, double ui){
  int m,n;
  VecGetOwnershipRange(u, &m, &n);
  if (i >= m && i < n){
    VecSetValues(u, 1, &i, &ui, ADD_VALUES); 
  }
}

void vecSetInsert(Vec u, int i, double ui){
  int m,n;
  VecGetOwnershipRange(u, &m, &n);
  if (i >= m && i < n){
    VecSetValues(u, 1, &i, &ui, INSERT_VALUES);
  }
}

void saveMat(Mat A, string filename){
  int world_rank; 
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank == 0) cout << "PETSc: Writing: " << filename << endl;
  PetscErrorCode code;

  if (get_filename_extension(filename) == "bin"){
    PetscViewer binaryViewer;
    code = PetscViewerBinaryOpen(MPI_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &binaryViewer); CHKERR(code);
    code = PetscViewerPushFormat(binaryViewer,PETSC_VIEWER_NATIVE); CHKERR(code);
    code = MatView(A, binaryViewer); CHKERR(code);
    code = PetscViewerDestroy(&binaryViewer); CHKERR(code);

  } else if (get_filename_extension(filename) == "m"){
    PetscViewer matlabViewer;
    code = PetscViewerCreate(MPI_COMM_WORLD,&matlabViewer); CHKERR(code);
    code = PetscViewerASCIIOpen(MPI_COMM_WORLD,filename.c_str(),&matlabViewer); CHKERR(code);
    code = PetscViewerPushFormat(matlabViewer,PETSC_VIEWER_ASCII_MATLAB); CHKERR(code);
    code = MatView(A, matlabViewer); CHKERR(code);
    code = PetscViewerDestroy(&matlabViewer); CHKERR(code);
  } else if (get_filename_extension(filename) == "txt"){
    PetscViewer matlabViewer;
    code = PetscViewerCreate(MPI_COMM_WORLD,&matlabViewer); CHKERR(code);
    code = PetscViewerASCIIOpen(MPI_COMM_WORLD,filename.c_str(),&matlabViewer); CHKERR(code);
    code = PetscViewerPushFormat(matlabViewer,PETSC_VIEWER_ASCII_MATLAB); CHKERR(code);
    code = MatView(A, matlabViewer); CHKERR(code);
    code = PetscViewerDestroy(&matlabViewer); CHKERR(code);
    ifstream txtFile(filename);
    string line;
    for (int i = 0; i < 6; i++){
      getline(txtFile, line);
    }
    ifstream txtFileOut(filename + ".txt");
    while (getline(txtFileOut, line)){
      cout << line << endl;
    }
    txtFileOut.close();
  }
}

void loadMat(Mat M, string filename){
  int world_rank; 
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank == 0) cout <<"PETSc: Reading: " << filename << endl; check_existence(filename);
  if (get_filename_extension(filename) == "bin"){
    PetscViewer binaryViewer;
    PetscErrorCode code;
    code = PetscViewerBinaryOpen(MPI_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &binaryViewer); CHKERR(code);
    code = PetscViewerPushFormat(binaryViewer,PETSC_VIEWER_NATIVE); CHKERR(code);
    code = MatLoad(M, binaryViewer); CHKERR(code); 
    code = PetscViewerDestroy(&binaryViewer); CHKERR(code);
  } else if (get_filename_extension(filename) == "txt"){
    vector<vector<double>> M_stl; 
    M_stl = importdata(filename);  
    for (int i = 0; i < M_stl.size(); i++){
      for (int j = 0; j < M_stl[0].size(); j++){
        MatSetValue(M, i, j, M_stl[i][j], INSERT_VALUES);
      }
    }
    MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  }
}


void loadVec(vector<double> & u, string filename){
  int world_rank; 
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank == 0) cout << "PETSc: Reading: " << filename; check_existence(filename);

  if (get_filename_extension(filename) == "vct"){
    int vecSize = u.size();

    ifstream ensightFile((filename).c_str());

    /* jump header */
    string line;
    getline(ensightFile,line);

    int iterVec = 0;
    if (ensightFile.is_open()){
      while (getline(ensightFile, line)){ 
        if (line.size() != 72)
          if (world_rank == 0) cout << "\nInvalid file " << filename << endl;
        u[iterVec] = stod(line.substr( 0, 12)); iterVec++;
        u[iterVec] = stod(line.substr(12, 12)); iterVec++;
        u[iterVec] = stod(line.substr(24, 12)); iterVec++;
        if (iterVec >= vecSize)
          break;
        u[iterVec] = stod(line.substr(36, 12)); iterVec++;
        u[iterVec] = stod(line.substr(48, 12)); iterVec++;
        u[iterVec] = stod(line.substr(60, 12)); iterVec++;
        if (iterVec >= vecSize)
          break;
      }
    } else {
      if (world_rank == 0) cout << "\nCan't open " << filename << endl;
      exit(1);
    }
  } else if (get_filename_extension(filename) == "scl"){
    int vecSize = u.size();

    ifstream ensightFile((filename).c_str());

    /* jump header */
    string line;
    getline(ensightFile,line);

    int iterVec = 0;
    if (ensightFile.is_open()){
      while (getline(ensightFile, line)){ 
        u[iterVec] = stod(line.substr( 0, 12)); iterVec++;
        if (iterVec >= vecSize)
          break;
        u[iterVec] = stod(line.substr(12, 12)); iterVec++;
        if (iterVec >= vecSize)
          break;
        u[iterVec] = stod(line.substr(24, 12)); iterVec++;
        if (iterVec >= vecSize)
          break;
        u[iterVec] = stod(line.substr(36, 12)); iterVec++;
        if (iterVec >= vecSize)
          break;
        u[iterVec] = stod(line.substr(48, 12)); iterVec++;
        if (iterVec >= vecSize)
          break;
        u[iterVec] = stod(line.substr(60, 12)); iterVec++;
        if (iterVec >= vecSize)
          break;
      }
    } else {
      if (world_rank == 0) cout << "Can't open " << filename << endl;
      exit(1);
    }
  } else {
    if (world_rank == 0) cout << "Filename extension not recognized." << endl;
    exit(1);
  }
}

void loadVec(Vec u, string filename){
  int world_rank; 
  PetscErrorCode code;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank == 0) cout << "PETSc: Reading: " << filename; check_existence(filename);
  if (get_filename_extension(filename) == "bin"){
    PetscViewer binaryViewer;
    PetscErrorCode code;
    code = PetscViewerBinaryOpen(MPI_COMM_WORLD, filename.c_str(), FILE_MODE_READ, &binaryViewer); CHKERR(code);
    code = PetscViewerPushFormat(binaryViewer,PETSC_VIEWER_NATIVE); CHKERR(code); 
    code = VecLoad(u, binaryViewer); CHKERR(code);
    code = PetscViewerDestroy(&binaryViewer); CHKERR(code);

  } else if (get_filename_extension(filename) == "vct"){
    int vecSize;
    code = VecGetSize(u, &vecSize); CHKERR(code);
    double u_array[vecSize];

    ifstream ensightFile((filename).c_str());

    /* jump header */
    string line;
    getline(ensightFile,line);

    int m,n;
    code = VecGetOwnershipRange(u, &m, &n); CHKERR(code);

    int iterVec = 0;
    if (ensightFile.is_open()){
      while (getline(ensightFile, line)){ 
        u_array[iterVec] = stod(line.substr( 0, 12)); iterVec++;
        u_array[iterVec] = stod(line.substr(12, 12)); iterVec++;
        u_array[iterVec] = stod(line.substr(24, 12)); iterVec++;
        if (iterVec >= vecSize)
          break;
        u_array[iterVec] = stod(line.substr(36, 12)); iterVec++;
        u_array[iterVec] = stod(line.substr(48, 12)); iterVec++;
        u_array[iterVec] = stod(line.substr(60, 12)); iterVec++;
        if (iterVec >= vecSize)
          break;
      }

      int * indexes = new int[n - m];
      double * u_arrayLocal = new double[n - m];
      int counter = 0;

      for (int i = m; i < n; i++){
        indexes[counter] = i;
        u_arrayLocal[counter] = u_array[i];
        counter++;
      }

      code = VecSetValues(u, n - m, indexes, u_arrayLocal, INSERT_VALUES); CHKERR(code);
      VecAssemblyBegin(u);
      VecAssemblyEnd(u);
    } else {
      if (world_rank == 0) cout << "\nCan't open " << filename << endl;
      exit(1);
    }
  } else if (get_filename_extension(filename) == "scl"){
    int vecSize;
    VecGetSize(u, &vecSize); 
    double u_array[vecSize];

    ifstream ensightFile((filename).c_str());

    /* jump header */
    string line;
    getline(ensightFile,line);

    int iterVec = 0;
    if (ensightFile.is_open()){
      while (getline(ensightFile, line)){ 
        u_array[iterVec] = stod(line.substr( 0, 12)); iterVec++;
        if (iterVec >= vecSize)
          break;
        u_array[iterVec] = stod(line.substr(12, 12)); iterVec++;
        if (iterVec >= vecSize)
          break;
        u_array[iterVec] = stod(line.substr(24, 12)); iterVec++;
        if (iterVec >= vecSize)
          break;
        u_array[iterVec] = stod(line.substr(36, 12)); iterVec++;
        if (iterVec >= vecSize)
          break;
        u_array[iterVec] = stod(line.substr(48, 12)); iterVec++;
        if (iterVec >= vecSize)
          break;
        u_array[iterVec] = stod(line.substr(60, 12)); iterVec++;
        if (iterVec >= vecSize)
          break;
      }

      int indexes[vecSize];
      for (int i = 0; i < vecSize; i++){
        indexes[i] = i;
      }
      VecSetValues(u, vecSize, indexes, u_array, INSERT_VALUES);
      VecAssemblyBegin(u);
      VecAssemblyEnd(u);
    } else {
      if (world_rank == 0) cout << "PETSc: Can't open " << filename << endl;
      exit(1);
    }
  } else {
    if (world_rank == 0) cout << "PETSc: Filename extension not recognized." << endl;
    exit(1);
  }
  double vecNorm;
  VecNorm(u, NORM_2, &vecNorm);
  if (world_rank == 0) cout << "\tNorm_l2: " << vecNorm << endl;
}

void saveVec(Vec u, string filename){
  int world_rank; 
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank == 0) cout << "PETSc: Writing: " << filename << endl;
  PetscViewer binaryViewer;
  PetscErrorCode code;
  code = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(), FILE_MODE_WRITE, &binaryViewer); CHKERR(code);
  code = PetscViewerPushFormat(binaryViewer,PETSC_VIEWER_NATIVE); CHKERR(code);
  code = VecView(u, binaryViewer); CHKERR(code);
  code = PetscViewerDestroy(&binaryViewer); CHKERR(code);
}

void vec(Vec & u, int sizeVec){
  PetscErrorCode code;
  code = VecCreate(PETSC_COMM_WORLD, &u); CHKERR(code);
  code = VecSetSizes(u, PETSC_DECIDE, sizeVec); CHKERR(code);
  code = VecSetFromOptions(u); CHKERR(code);
  code = VecZeroEntries(u); CHKERR(code);
}

vector<Vec> vec(vector<Vec> u){
  PetscErrorCode code;
  int sizeVec;
  vector<Vec> v;
  v.resize(u.size());
  for (int i = 0 ; i < u.size(); i++){
    code = VecGetSize(u[i], &sizeVec); CHKERR(code);
    code = VecCreate(PETSC_COMM_WORLD, &v[i]); CHKERR(code);
    code = VecSetSizes(v[i], PETSC_DECIDE, sizeVec); CHKERR(code);
    code = VecSetFromOptions(v[i]); CHKERR(code);
    code = VecCopy(u[i], v[i]); CHKERR(code);
  }
  return v;
}

Vec vec(Vec u){
  PetscErrorCode code;
  int sizeVec;
  Vec v;
  code = VecGetSize(u, &sizeVec); CHKERR(code);
  code = VecCreate(PETSC_COMM_WORLD, &v); CHKERR(code);
  code = VecSetSizes(v, PETSC_DECIDE, sizeVec); CHKERR(code);
  code = VecSetFromOptions(v); CHKERR(code);
  code = VecCopy(u, v); CHKERR(code);
  return v;
}

Vec vec(int sizeVec, string path_vec){
  Vec u;
  vec(u, sizeVec);
  loadVec(u, path_vec);
  return u;
}

Vec vec(int sizeVec){
  Vec u;
  vec(u, sizeVec);
  return u;
}

void mat(Mat & A, int nbRows, int nbCols, string type){
  assert(type == "sparse" or type == "dense");
  PetscErrorCode code;
  if (type == "sparse"){
    MatCreate(PETSC_COMM_WORLD, &A); 
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, nbRows, nbCols);
  } else if (type == "dense"){
    MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nbRows, nbCols, NULL, &A);
  }
  code = MatSetFromOptions(A); CHKERR(code); 
  code = MatSetUp(A); CHKERR(code); 
  code = MatZeroEntries(A); CHKERR(code);
}

Mat mat(int nbRows, int nbCols, string type){
  Mat A;
  assert(type == "sparse" or type == "dense");
  if (type == "sparse"){
    MatCreate(PETSC_COMM_WORLD, &A); 
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, nbRows, nbCols);
  } else if (type == "dense"){
    MatCreateDense(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, nbRows, nbCols, NULL, &A);
  } else {
    exit(1);
  }
  MatSetFromOptions(A); 
  MatSetUp(A); 
  MatZeroEntries(A);
  return A;
}

Mat mat(Mat A){
  PetscErrorCode code;
  int m, n;
  MatGetSize(A, &m, &n);
  Mat B = mat(m, n);
  code = MatDuplicate(A, MAT_COPY_VALUES, &B); CHKERR(code);
  MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
  return B;
}


void configureKSP(KSP & ksp, Mat & A, string solver, string preconditioner, bool monitor, bool use_solution_as_guess_KSP, bool reuse_preconditioner, double ksp_tolerance){

  int world_rank; 
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  PetscErrorCode code;
  code = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERR(code); 
  code = KSPSetOperators(ksp, A, A); CHKERR(code); 
  if (solver == "gmres"){
    code = KSPSetType(ksp, KSPGMRES); CHKERR(code); 
  } else if (solver == "preonly"){
    code = KSPSetType(ksp, KSPPREONLY); CHKERR(code); 
  } else if (solver == "cg"){
    code = KSPSetType(ksp, KSPCG); CHKERR(code); 
  } else {
//    if (world_rank == 0) cout << "Solver method " << solver << " not valid." << endl;
//    exit(1);
  }
  PC pc;
  code = KSPGetPC(ksp, &pc); CHKERR(code); 
  if (preconditioner == "asm"){
    code = PCSetType(pc, PCASM); CHKERR(code); 
  } else if (preconditioner == "gasm"){
    code = PCSetType(pc, PCGASM); CHKERR(code); 
  } else if (preconditioner == "lu"){
    code = PCSetType(pc, PCLU); CHKERR(code); 
  } else if (preconditioner == "ilu"){
    code = PCSetType(pc, PCILU); CHKERR(code); 
  } else if (preconditioner == "jacobi"){
    code = PCSetType(pc, PCJACOBI); CHKERR(code); 
  } else if (preconditioner == "none"){
    code = PCSetType(pc, PCNONE); CHKERR(code); 
  } else {
//    if (world_rank == 0) cout << "Preconditioner " + preconditioner + " not valid." << endl;
//    exit(1);
  }
  if (monitor && solver != "preonly"){
    code = KSPMonitorSet(ksp, krylovMonitor, NULL, NULL); CHKERR(code);
  }
  if (use_solution_as_guess_KSP && solver != "preonly"){
    code = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERR(code);
  }
  if (reuse_preconditioner){
    code = KSPSetReusePreconditioner(ksp, PETSC_TRUE); CHKERR(code);
  }
  code = KSPSetTolerances(ksp, ksp_tolerance, 1e-50, 20000, 20000);
  code = KSPSetFromOptions(ksp); CHKERR(code); 
}

void configureKSP(KSP & ksp, Parameters par){

  int world_rank; 
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  PetscErrorCode code;
  code = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERR(code); 
  if (par.solver() == "gmres"){
cout << "ASDASD" << endl;
    code = KSPSetType(ksp, KSPGMRES); CHKERR(code); 
    code = KSPGMRESSetRestart(ksp, par.ksp_restartGMRESiterations()); CHKERR(code);
  } else if (par.solver() == "fgmres"){
    code = KSPSetType(ksp, KSPFGMRES); CHKERR(code); 
    code = KSPGMRESSetRestart(ksp, par.ksp_restartGMRESiterations()); CHKERR(code);
  } else if (par.solver() == "preonly"){
    code = KSPSetType(ksp, KSPPREONLY); CHKERR(code); 
  } else if (par.solver() == "cg"){
    code = KSPSetType(ksp, KSPCG); CHKERR(code); 
  } else {
    if (world_rank == 0) cout << "Solver method " << par.solver() << " not valid." << endl;
    exit(1);
  }
  PC pc;
  code = KSPGetPC(ksp, &pc); CHKERR(code); 
  if (par.preconditioner() == "asm"){
cout << "PREPRE" << endl;
    code = PCSetType(pc, PCASM); CHKERR(code); 
  } else if (par.preconditioner() == "gasm"){
    code = PCSetType(pc, PCGASM); CHKERR(code); 
  } else if (par.preconditioner() == "amg"){
    code = PCSetType(pc, PCGAMG); CHKERR(code); 
  } else if (par.preconditioner() == "lu"){
    code = PCSetType(pc, PCLU); CHKERR(code); 
  } else if (par.preconditioner() == "cholesky"){
    code = PCSetType(pc, PCCHOLESKY); CHKERR(code); 
  } else if (par.preconditioner() == "ilu"){
    code = PCSetType(pc, PCILU); CHKERR(code); 
  } else if (par.preconditioner() == "jacobi"){
    code = PCSetType(pc, PCJACOBI); CHKERR(code); 
  } else if (par.preconditioner() == "vanka"){
    code = PCSetType(pc, PCPATCH); CHKERR(code); 
  } else if (par.preconditioner() == "none"){
    code = PCSetType(pc, PCNONE); CHKERR(code); 
  } else {
    if (world_rank == 0) cout << "Preconditioner " + par.preconditioner() + " not valid." << endl;
    exit(1);
  }
  if (par.monitorKSP()){
    code = KSPMonitorSet(ksp, krylovMonitor, NULL, NULL); CHKERR(code);

  }
  if (par.use_solution_as_guess_KSP() && par.solver() != "preonly"){
    code = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERR(code);
  }

  if (par.solver() == "gmres" and par.useModifiedGramSchmidt() == true){
    KSPGMRESSetOrthogonalization(ksp, KSPGMRESModifiedGramSchmidtOrthogonalization);
  }
  if (par.reuse_preconditioner()){
    code = KSPSetReusePreconditioner(ksp, PETSC_TRUE); CHKERR(code);
  }
  code = KSPSetTolerances(ksp, par.ksp_tolerance(), par.ksp_tolerance_absolute(), PETSC_DEFAULT, par.ksp_max_iterations());
  code = KSPSetFromOptions(ksp); CHKERR(code); 
}

void configureKSP(KSP & ksp, string solver, string preconditioner, bool monitor, bool use_solution_as_guess_KSP, bool reuse_preconditioner, double ksp_tolerance){

//  cout << "this configureKSP function is bugged " << endl;
//  exit(1);
  int world_rank; 
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  PetscErrorCode code;
  code = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERR(code); 
  if (solver == "gmres"){
    code = KSPSetType(ksp, KSPGMRES); CHKERR(code); 
  } else if (solver == "preonly"){
    code = KSPSetType(ksp, KSPPREONLY); CHKERR(code); 
  } else if (solver == "cg"){
    code = KSPSetType(ksp, KSPCG); CHKERR(code); 
  } else {
    if (world_rank == 0) cout << "Solver method " << solver << " not valid." << endl;
    exit(1);
  }
  PC pc;
  code = KSPGetPC(ksp, &pc); CHKERR(code); 
  if (preconditioner == "asm"){
    code = PCSetType(pc, PCASM); CHKERR(code); 
  } else if (preconditioner == "lu"){
    code = PCSetType(pc, PCLU); CHKERR(code); 
  } else if (preconditioner == "cholesky"){
    code = PCSetType(pc, PCCHOLESKY); CHKERR(code); 
  } else if (preconditioner == "ilu"){
    code = PCSetType(pc, PCILU); CHKERR(code); 
  } else if (preconditioner == "jacobi"){
    code = PCSetType(pc, PCJACOBI); CHKERR(code); 
  } else if (preconditioner == "vanka"){
    code = PCSetType(pc, PCPATCH); CHKERR(code); 
  } else if (preconditioner == "none"){
    code = PCSetType(pc, PCNONE); CHKERR(code); 
  } else {
    if (world_rank == 0) cout << "Preconditioner " + preconditioner + " not valid." << endl;
    exit(1);
  }
  if (monitor){
    code = KSPMonitorSet(ksp, krylovMonitor, NULL, NULL); CHKERR(code);
  }
  if (use_solution_as_guess_KSP){
    code = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE); CHKERR(code);
  }
  if (reuse_preconditioner){
    code = KSPSetReusePreconditioner(ksp, PETSC_TRUE); CHKERR(code);
  }
  code = KSPSetTolerances(ksp, ksp_tolerance, 1e-50, 20000, 20000); CHKERR(code);
  code = KSPSetFromOptions(ksp); CHKERR(code); 
}

void configureEPS(Parameters par, EPS eps, Mat A){
  PetscInt nconv;
  EPSCreate(PETSC_COMM_WORLD, &eps);
  EPSSetOperators(eps, A, NULL);
  EPSSetProblemType(eps, EPS_GHEP);
  EPSSetType(eps, EPSKRYLOVSCHUR);
  EPSSetDimensions(eps, par.nbModes(), PETSC_DECIDE, PETSC_DECIDE);
  EPSSetFromOptions(eps);
}

Vec zeros(int nbZeros){
  PetscErrorCode code;
  Vec u;
  vec(u, nbZeros);
  code = VecZeroEntries(u); CHKERR(code);
  return u;
}

void zeros(Vec & u, int nbZeros){
  vec(u, nbZeros);
  VecZeroEntries(u);
}

KSP configureKSP(Mat & A, string solver, string preconditioner, bool monitor){
  int world_rank; 
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank == 0) cout << "PETSc: Configuring Krylov solver" << endl;
  KSP ksp;
  configureKSP(ksp, A, solver, preconditioner, monitor);
  
  return ksp;
}

PetscErrorCode krylovMonitor(KSP ksp,int n,PetscReal rnorm,void *dummy){
  PetscErrorCode code;
  code = PetscPrintf(PETSC_COMM_WORLD,"    KSP: iteration %D, residual norm:  %14.12e \r",n,rnorm);CHKERR(code);
  return 0;
}

/* function kept for compatibility, use krylovMonitor */
PetscErrorCode MyKSPMonitor(KSP ksp,int n,PetscReal rnorm,void *dummy){
  PetscErrorCode ierr;
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\titeration %D KSP Residual norm %14.12e \n",n,rnorm);CHKERRQ(ierr);
  return 0;
}

vector<double> getSingularValues(Mat A, int n){
  vector<double> sValues;
  sValues.resize(n);
  SVD svd;           
  SVDCreate(PETSC_COMM_WORLD, &svd);  
  SVDSetOperator(svd,A);  
  SVDSetDimensions(svd, n, PETSC_DEFAULT, PETSC_DEFAULT);  
  SVDSetFromOptions(svd);  
  SVDSolve(svd);  
  for (int i = 0; i < n; i++){
    SVDGetSingularTriplet(svd, i, &sValues[i], NULL, NULL); 
  }
  SVDDestroy(&svd); 
  return sValues;
}

void CHKERR(PetscErrorCode code){
  if (code != 0)
    exit(1);
}

double norm(Vec u){
  PetscErrorCode code;
  double norm;
  code = VecNorm(u, NORM_2, &norm); CHKERR(code);
  return norm;
}

double norm(Mat m){
  PetscErrorCode code;
  PetscScalar norm;
  code = MatNorm(m, NORM_FROBENIUS, &norm); CHKERR(code);
  return norm;
}

double dot(Vec u, Vec v){
  PetscErrorCode code;
  double dotp;
  code = VecDot(u, v, &dotp); CHKERR(code);
  return dotp;
}

ostream& operator<<(ostream& out, Vec u){
  int world_rank; 
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  PetscErrorCode code;
  int sizeVec;
  code = VecGetSize(u, &sizeVec); CHKERR(code);
  double * u_array = new double[sizeVec];
  int indexes[sizeVec];
  for (int i = 0; i < sizeVec; i++){
    indexes[i] = i;
  }
  code = VecGetValues(u, sizeVec, indexes, u_array); CHKERR(code); 
  if (world_rank == 0) cout << "[ ";
  for (int i = 0; i < sizeVec; i++){
    out << u_array[i] << " ";
  }
  if (world_rank == 0) cout << " ]";
  return out;
}

//Mat operator*(Mat A, Mat B){
//  int m, n;
//  MatGetSize(A, &m, NULL);
//  MatGetSize(B, NULL, &n);
//  Mat P = mat(m, n);
//  MatMatMult(A, B, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &P);
//  return P;
//}

ostream& operator<<(ostream& out, Mat m){
  PetscErrorCode code;
  int nbRows,nbCols;
  code = MatGetSize(m, &nbRows, &nbCols); CHKERR(code);
  assert(nbRows < 50 and nbCols < 50); /* This function only makes sense for non memory-consuming matrices */
  vector<double> m_ij(nbCols * nbRows);
  MatGetValues(m, nbRows, &range(nbRows)[0], nbCols, &range(nbCols)[0], &m_ij[0]);
  for (int i : range(nbRows)){
    out << "[ ";
    for (int j : range(nbCols)){
      out << scientific << m_ij[nbRows*i + j] << " ";    
    }
    out << " ]\n";
  }
  return out;
}

double ip(Vec u, Mat M, Vec v){
  int size;
  VecGetSize(u, &size);
  PetscReal ip_val;
  Vec mtv;
  vec(mtv, size);
  MatMult(M, v, mtv);
  VecDot(u, mtv, &ip_val);
  VecDestroy(&mtv);
  return ip_val;
}

void blockMatrix(Mat & A, vector<int> index){
  PetscErrorCode code;
  code = MatZeroRows(A, index.size(), &index[0], 1.0, NULL, NULL); CHKERR(code);
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERR(code);
}

Mat outer(Vec u){
  int N;
  VecGetSize(u, &N);
  Mat C = mat(N, N, "dense");
  double u_arr[N];
  VecGetValues(u, N, &range(N)[0], &u_arr[0]);
  PetscErrorCode code;
  for (int i=0;i<N;i++){
    for (int j=0;j<N;j++){
      code = MatSetValue(C,i,j,u_arr[i]*u_arr[j],INSERT_VALUES); CHKERR(code);
    }
  }
  code = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY); CHKERR(code);
  return C;
}

Mat diag(vector<double> d){
  Mat A = mat(d.size(), d.size());
  PetscErrorCode code;
  code = MatDiagonalSet(A, petsc(d), INSERT_VALUES); CHKERR(code);
  return A;
}

/* Place the components of basis in the columns of a matrix */
Mat buildProjector(const vector<Vec> & basis){
  PetscErrorCode code;
  int nCols = basis.size();
  int nRows;
  VecGetSize(basis[0], &nRows);
  Mat A = mat(nRows, nCols, "dense");

  int m,n;
  MatGetOwnershipRange(A, &m, &n); 
  int nbDofsLocal = n - m;
  for (int i=0;i<nCols;i++){
    double col[nbDofsLocal];
    code = VecGetValues(basis[i], nbDofsLocal, &range(m, n)[0], col); CHKERR(code);
    code = MatSetValues(A, nbDofsLocal, &range(m,n)[0], 1, &i, col, INSERT_VALUES); CHKERR(code);
  }
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  return A;
}

Vec petsc(vector<double> u){
  Vec v = vec(u.size());
  PetscErrorCode code;
  code = VecSetValues(v, u.size(), &range(u.size())[0], &u[0], INSERT_VALUES); CHKERR(code);
  code = VecAssemblyBegin(v); CHKERR(code);
  code = VecAssemblyEnd(v); CHKERR(code);
  return v;
}

Mat petsc(vector<vector<double>> m){
  Mat A = mat(m.size(), m[0].size(), "dense");
  PetscErrorCode code;
  for (int i = 0; i < m.size(); i++){
    for (int j = 0; j < m[i].size(); j++){
      code = MatSetValue(A, i, j, m[i][j], INSERT_VALUES); CHKERR(code);
    }
  }
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  return A;
}

vector<double> stl(Vec u){
  PetscErrorCode code;
  int N;
  Vec u_seq = getSequential(u);
  code = VecGetSize(u_seq, &N); CHKERR(code);
  vector<double> v(N);
  code = VecGetValues(u_seq, N, &range(N)[0], &v[0]); CHKERR(code);
  code = VecDestroy(&u_seq); CHKERR(code);
  return v;
}

Vec ones(int n){
  Vec one = vec(n);
  for (int i=0; i<n; i++){
    double uno = 1.0;
    VecSetValues(one, 1, &i, &uno, INSERT_VALUES);
  }
  VecAssemblyBegin(one);
  VecAssemblyEnd(one);
  return one;
}

Mat eye(int n, string type){
  Mat I = mat(n, n, type);
  Vec D = ones(n);
  MatDiagonalSet(I, D, INSERT_VALUES);
  return I;
}

Mat transpose(Mat A){
  PetscErrorCode code;
  int n, m;
  MatGetSize(A, &m, &n);
  Mat At = mat(n, m);
  code = MatTranspose(A, MAT_INITIAL_MATRIX, &At); CHKERR(code);
  return At;
}

Mat inverse(Mat A){
  PetscErrorCode code;
  int m,n;
  MatGetSize(A, &m, &n);

  Mat U = mat(m,m, "dense");
  Mat V = mat(m,m, "dense");
  Mat S = mat(m,m, "dense");
  SVD svd;           
  SVDCreate(PETSC_COMM_WORLD, &svd);  
  SVDSetOperator(svd, A);  
  SVDSetDimensions(svd, m, PETSC_DEFAULT, PETSC_DEFAULT);  
  SVDSetFromOptions(svd);  
  SVDSolve(svd);  
  for (int i = 0; i < m; i++){
    double svalue;
    Vec u = vec(m);
    Vec v = vec(m);
    SVDGetSingularTriplet(svd, i, &svalue, u, v); 
    MatSetValues(V, m, &range(n)[0], 1, &i, &stl(v)[0], INSERT_VALUES);
    MatSetValues(U, m, &range(n)[0], 1, &i, &stl(u)[0], INSERT_VALUES);
    MatSetValue(S, i, i, 1.0/(svalue), INSERT_VALUES);
  }
  MatAssemblyBegin(V, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(V, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(U, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(U, MAT_FINAL_ASSEMBLY);

  Mat Ut = mat(m, m);  
  code = MatTranspose(U, MAT_INITIAL_MATRIX, &Ut); CHKERR(code);

  Mat B = mat(m,m,"dense");
  MatMatMult(V, S, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &B);

  Mat C = mat(m,m,"dense");
  MatMatMult(B, Ut, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C);

  SVDDestroy(&svd); 
  return C;
}

vector<vector<double>> stl(Mat A){
  int m, n;
  MatGetSize(A, &m, &n);
  vector<vector<double>> A_stl(m);
  for (int i = 0; i < m; i++){
    A_stl.resize(n);
    for (int j = 0; j < n; j++){
      vector<double> value(1);
      MatGetValues(A, 1, &i, 1, &j, &value[0]);
      A_stl[i][j] = value[0];
    }
  }
  return A_stl;
}

/* Compute ATA = V * invS * Ut */
Mat bestATAinverse(Mat A){

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank == 0) cout << "PETSc: Computing matrix pseudo-inverse." << endl;

  PetscErrorCode code;
  int m,n, nconv;
  MatGetSize(A, &m, &n);

  SVD svd;           
  SVDCreate(PETSC_COMM_WORLD, &svd);  
  SVDSetOperator(svd, A);  
  SVDSetDimensions(svd, m, PETSC_DEFAULT, PETSC_DEFAULT);  
  SVDSetFromOptions(svd);  
  SVDSolve(svd);  
  SVDGetConverged(svd, &nconv);
  vector<double> csv = getSingularValues(A, nconv); // converged singular values
  vector<double> nzsv; // non-zero singular values
  for (double sv : csv){
    if (sv > 0.000001){ // 10^{-6}
      nzsv.push_back(sv); 
    }
  }
  int nnzsv = nzsv.size();

  if (world_rank == 0) cout << "PETSc: Number of non-zero singular values: " << nnzsv << "." << endl;
  if (world_rank == 0 && nnzsv < n) cout << "PETSc: Using SVD regularization." << endl;
  if (world_rank == 0) cout << "PETSc: Singular values: " << nzsv << "." << endl;
  
  Mat U = mat(m, nnzsv, "dense");
  Mat V = mat(n, nnzsv, "dense");
  Mat S = mat(nnzsv, nnzsv, "dense");
  for (int i = 0; i < nnzsv; i++){
    Vec u_vec = vec(m);
    Vec v_vec = vec(n);
    SVDGetSingularTriplet(svd, i, NULL, u_vec, v_vec); 
    MatSetValues(U, m, &range(m)[0], 1, &i, &stl(u_vec)[0], INSERT_VALUES);
    MatSetValues(V, n, &range(n)[0], 1, &i, &stl(v_vec)[0], INSERT_VALUES);
    MatSetValue(S, i, i, 1.0/nzsv[i], INSERT_VALUES);
  }
  MatAssemblyBegin(V, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(V, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(U, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(U, MAT_FINAL_ASSEMBLY);

  Mat Ut = mat(nnzsv, m);  
  code = MatTranspose(U, MAT_INITIAL_MATRIX, &Ut); CHKERR(code);

  Mat B = mat(n,nnzsv,"dense");
  MatMatMult(V, S, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &B);

  Mat C = mat(n,m,"dense");
  MatMatMult(B, Ut, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C);

  SVDDestroy(&svd); 
  return C;
}

/* Compute AAT = U * invS * Vt */
Mat bestAATinverse(Mat A){

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank == 0) cout << "PETSc: Computing matrix pseudo-inverse." << endl;

  PetscErrorCode code;
  int m,n, nconv;
  MatGetSize(A, &m, &n);

  SVD svd;           
  SVDCreate(PETSC_COMM_WORLD, &svd);  
  SVDSetOperator(svd, A);  
  SVDSetDimensions(svd, m, PETSC_DEFAULT, PETSC_DEFAULT);  
  SVDSetFromOptions(svd);  
  SVDSolve(svd);  
  SVDGetConverged(svd, &nconv);
  vector<double> csv = getSingularValues(A, nconv); // converged singular values
  vector<double> nzsv; // non-zero singular values
  for (double sv : csv){
    if (sv > 0.000001){ // 10^{-6}
      nzsv.push_back(sv); 
    }
  }
  int nnzsv = nzsv.size();

  if (world_rank == 0) cout << "PETSc: Number of non-zero singular values: " << nnzsv << "." << endl;
  if (world_rank == 0 && nnzsv < m) cout << "PETSc: Using SVD regularization." << endl;

  Mat U = mat(m, nnzsv, "dense");
  Mat V = mat(n, nnzsv, "dense");
  Mat S = mat(nnzsv, nnzsv, "dense");
  for (int i = 0; i < nnzsv; i++){
    Vec u_vec = vec(m);
    Vec v_vec = vec(n);
    SVDGetSingularTriplet(svd, i, NULL, u_vec, v_vec); 
    MatSetValues(U, m, &range(m)[0], 1, &i, &stl(u_vec)[0], INSERT_VALUES);
    MatSetValues(V, n, &range(n)[0], 1, &i, &stl(v_vec)[0], INSERT_VALUES);
    MatSetValue(S, i, i, 1.0/nzsv[i], INSERT_VALUES);
  }
  MatAssemblyBegin(V, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(V, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY);
  MatAssemblyBegin(U, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(U, MAT_FINAL_ASSEMBLY);

  Mat Vt = mat(nnzsv, n, "dense");  
  code = MatTranspose(V, MAT_INITIAL_MATRIX, &Vt); CHKERR(code);

  Mat B = mat(m, nnzsv,"dense");
  MatMatMult(U, S, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &B);

  Mat C = mat(m,n,"dense");
  MatMatMult(B, Vt, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &C);
  SVDDestroy(&svd); 
  return C;

}

Vec getSequential(Vec u){
  int u_size;
  PetscErrorCode code;
  VecGetSize(u, &u_size);

  /* scatter to sequential vector */
  Vec u_seq;
  code = VecCreateSeq(PETSC_COMM_SELF, u_size, &u_seq); CHKERR(code);
  code = VecSetFromOptions(u_seq); CHKERR(code); 
  VecScatter inctx;
  code = VecScatterCreateToAll(u, &inctx, NULL); CHKERR(code);
  code = VecScatterBegin(inctx, u, u_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERR(code);

  code = VecScatterEnd(inctx, u, u_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERR(code);
  code = VecScatterDestroy(&inctx); CHKERR(code);
  return u_seq;
}

//vector<Vec> split(Vec u, vector<int> nbDofsPerNode, int nbVertices){
//  PetscErrorCode code;
//  vector<Vec> split_vector(nbDofsPerNode.size());
//  for (int i = 0; i < nbDofsPerNode.size(); i++){
//    int offset = 0;
//    for (int j = 0; j < i; j++){
//      offset = offset + nbDofsPerNode[i]*nbVertices;
//    }
//    int nbDofs = nbDofsPerNode[i]*nbVertices;
//    double * sub_vector = new double [nbDofs];
//
//    code = VecGetValues(u, nbDofs, &range(offset, offset + nbDofs)[0], sub_vector); CHKERR(code);
//    vec(split_vector[i], nbDofs);
//    code = VecSetValues(split_vector[i], nbDofs, &range(nbDofs)[0], sub_vector, INSERT_VALUES); CHKERR(code);
//    code = VecAssemblyBegin(split_vector[i]); CHKERR(code);
//    code = VecAssemblyEnd(split_vector[i]); CHKERR(code);
//  }
//  return split_vector;
//}

//void split(vector<Vec> split_vector, Vec u, vector<int> nbDofsPerNode, int nbVertices){
//  PetscErrorCode code;
//  split_vector.resize(nbDofsPerNode.size());
//  for (int i = 0; i < nbDofsPerNode.size(); i++){
//    int offset = 0;
//    for (int j = 0; j < i; j++){
//      offset = offset + nbDofsPerNode[i]*nbVertices;
//    }
//    int nbDofs = nbDofsPerNode[i]*nbVertices;
//    double * sub_vector = new double [nbDofs];
//    code = VecGetValues(u, nbDofs, &range(offset, offset + nbDofs)[0], sub_vector); CHKERR(code);
//    vec(split_vector[i], nbDofs);
//    code = VecSetValues(split_vector[i], nbDofs, &range(nbDofs)[0], sub_vector, INSERT_VALUES); CHKERR(code);
//    code = VecAssemblyBegin(split_vector[i]); CHKERR(code);
//    code = VecAssemblyEnd(split_vector[i]); CHKERR(code);
//  }
//}

/* compute a^T M b */
double quadraticForm(Vec a, Mat M, Vec b){
  int m,n,k,l;
  VecGetSize(a, &m);
  VecGetSize(b, &n);
  MatGetSize(M, &k, &l);
  assert(m == k && n == l);

  double qf;
  Vec Mb = vec(n);
  MatMult(M, b, Mb);
  VecDot(Mb, a, &qf);
  return qf;
}
