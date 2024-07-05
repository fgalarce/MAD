/*=============================================================================
  This file is part of the code MAD
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2017 - 2021,
    
     Felipe Galarce at WIAS

  MDMA is free software; you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by the Free
  Software Foundation; either version 2.1 of the License, or (at your option)
  any later version.

  MDMA is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with MDMA. If not, see http://www.gnu.org/licenses/.
  =============================================================================*/

#include <mad.hpp>

int main(int argc, char *argv[]){

  assert(argc == 2);

  string data_file = argv[1];
  
  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  Parameters par(data_file);
  par.print();

  /* Initialize MDMA objects */
  IO io;
  io.initialize(par);

  Geometry geo;
  geo.initialize(par, io);

  InnerProduct ip;
  ip.initialize(par, geo);

  LinearAlgebra la;
  la.initialize(par, ip);

  Measures measures;
  measures.initialize(par, geo, la, io);

  Model Vn;
  Vn.initialize(par, geo, la, io);

  int m = measures.nbMeasures();
  int n = par.nbModes();

  int nbVertices = geo.nbVertices;
  int nbDofs = par.nbDofsPerNode()[0]*geo.nbVertices; 

//  if (par.rieszRepresentersFolder() == ""){
//    IO io_measures;
//    io_measures.initialize(par, "Observer");
//    for (int i = 0; i < m; i++){
//      io_measures.writeState(measures.rieszRepresenters(i), "observer", i);
//    }
//  }

//  IO io_model;
//  io_model.initialize(par, "Model");
//  for (int i = 0; i < n; i++){
//    io_model.writeState(Vn.basis(0.0)[i], "pod", i);
//  }

  if (world_rank == 0) cout << "PBDW: assembling cross-Gramian G(Wm, Vn)." << endl;
  vector<vector<double>> G_array(m);
  for (int i = 0; i < m; i++){ 
    G_array[i].resize(n);
    for (int j = 0; j < n; j++){
      G_array[i][j] = ip(Vn.basis(0.0)[j], measures.rieszRepresenters(i));
    }
    if (world_rank == 0) cout << "PBDW: row " << i << " norm: " << norm(G_array[i]) << endl;
  }

  /* Place Riesz representers in the columns of a matrix */
  Mat W = mat(m, nbDofs);
  Mat W_tr = mat(nbDofs, m);
  W = buildProjector(measures.rieszRepresenters());
  W_tr = transpose(W);

  Mat G = mat(m,n);
  Mat Gt = mat(n,m);
  int low, high;
  code = MatGetOwnershipRange(G, &low, &high); CHKERRQ(code);
  if (world_rank == 0) cout << m << " " << n << endl;
  for (int i = 0; i < m; i++){
    for (int j = 0; j < n; j++){
      if (i >= low && i < high){
        code = MatSetValue(G, i, j, G_array[i][j], INSERT_VALUES); CHKERRQ(code);
      }
    }
  }
  code = MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  MatTranspose(G, MAT_INITIAL_MATRIX, &Gt);

  Mat L = mat(n,n);
  MatMatMult(Gt, G, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &L);
  KSP ksp = configureKSP(L, "preonly", "lu");
  
  vector<double> svCrosGram = getSingularValues(L, n);
  exportData(par.dirResults() + "/svGramian.txt", svCrosGram);
  ofstream errorFile(par.dirResults() + "/error.txt");

  for (int iteration = par.start(); iteration < par.end(); iteration++ ){

    if (world_rank == 0) cout << "PBDW: - - - - - - - - - - - - - - - - - - - - - " << endl;
    if (world_rank == 0) cout << "PBDW: Reconstruction iteration: " << iteration<< endl;

    measures.computeSynthetic(iteration);
    Vec rhs = vec(n);
    MatMult(Gt, measures.measures(), rhs);
    Vec c = vec(n);
    KSPSolve(ksp, rhs, c);

    Vec u_star = zeros(nbDofs);
    vector<double> ci = stl(c);
    for (int i = 0; i < n; i++){
      VecAXPY(u_star, ci[i], Vn.basis(0)[i]);
    }

    Vec Measure = zeros(nbDofs);
    for (int i = 0; i < m; i++){
      VecAXPY(Measure, stl(measures.measures())[i], measures.rieszRepresenters()[i]);
    }

    if (world_rank == 0) cout << "PBDW: Computing corrector in W_m cap V_n^perp." << endl;
    Vec Wt_times_v = vec(m);
    Vec z = vec(nbDofs);
    code = MatMult(W, measures.measures(), z); CHKERRQ(code);
    Vec z1 = vec(nbDofs);
    if (par.innerProduct() != "l2"){
      Vec M_u_star = vec(nbDofs);
      code = MatMult(ip.matrix, u_star, M_u_star); CHKERRQ(code);
      code = MatMult(W_tr, M_u_star, Wt_times_v); CHKERRQ(code);
    } else {
      code = MatMult(W_tr, u_star, Wt_times_v); CHKERRQ(code);
    }
    code = MatMult(W, Wt_times_v, z1); CHKERRQ(code);
    code = VecAXPY(z, -1.0, z1); CHKERRQ(code);
    code = VecAXPY(u_star, 1.0, z); CHKERRQ(code);
  
    if (world_rank == 0) cout << "PBDW: Comparing with ground truth." << endl;
    Vec u = vec(nbDofs);
    io.loadVector(u, par.dirSyntheticField() + "." + wildcard(iteration) + ".vct");

    Vec error = zeros(nbDofs);
    VecAXPY(error, 1.0, u);
    VecAXPY(error, -1.0, u_star);

    double normU = sqrt(ip(u, u));
    double normError = sqrt(ip(error, error));

    if (par.innerProduct() == "H1"){
      double normGradU = sqrt(ip.only_grad(u, u));
      double normGradError = sqrt(ip.only_grad(error, error));
      errorFile << scientific << normError << " " << normU << " " << normGradError << " " << normGradU << endl; 
      if (world_rank == 0) cout << "PBDW: relative error gradient: " << normGradError / normGradU << endl;
    } else {
      errorFile << scientific << normError << " " << normU << endl; 
    }

    if (world_rank == 0) cout << "PBDW: norm error = " << normError << endl;
    if (world_rank == 0) cout << "PBDW: relative error = " << normError / normU << endl;
  
    io.writeState(error, "error", iteration);
    io.writeState(u_star, "u_star", iteration);
    io.writeState(Measure, "measures", iteration);
    io.writeState(u, "u", iteration);
    io.writeState(z, "update", iteration);
    VecDestroy(&Measure);
    VecDestroy(&u_star);
  }
  errorFile.close();
  measures.finalize();
  par.finalize();
  SlepcFinalize();
}
