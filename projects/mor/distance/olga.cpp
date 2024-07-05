/*=============================================================================
  This file is part of the code M.D.M.A. (or MDMA)
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2020,
    
     Felipe Galarce at INRIA

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
    loadVec(phi1[i], par.model() + "/mode." + wildcard(i) + ".vct");
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

  Mat M = mat(2*n, 2*n);
  cout << "MAIN: assembling M." << endl;
  for (int i=0; i<n; i++){
    for (int j=0; j<n; j++){
      code = MatSetValue(M, i, j, ip(phi0[i], phi0[j]), INSERT_VALUES); CHKERRQ(code);
      code = MatSetValue(M, n+i, n+j, ip(phi1[i], phi1[j]), INSERT_VALUES); CHKERRQ(code);
      code = MatSetValue(M, i, n+j, ip(phi0[i], phi1[j]), INSERT_VALUES); CHKERRQ(code);
      code = MatSetValue(M, n+i, j, ip(phi0[i], phi1[j]), INSERT_VALUES); CHKERRQ(code);
    }
  }
  code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

  Mat A = mat(2*n, 2*n);
  double Aij;
  for (int i=0; i<n; i++){
  cout << "MAIN: assembling A. Row: " << i << "."<< endl;
    for (int j=i; j<n; j++){
      Vec e_P_E_fi = zeros(nbDofs); /* e_i - P_F e_i */
      Vec e_P_E_fj = zeros(nbDofs); /* e_j - P_F e_j */
      Vec e_P_F_ei = zeros(nbDofs); /* e_i - P_F e_i */
      Vec e_P_F_ej = zeros(nbDofs); /* e_j - P_F e_j */
      for (int k=0; k<n; k++){
        code = VecAXPY(e_P_F_ei, -ip(phi0[i], phi1[k]), phi1[k]); CHKERRQ(code); 
        code = VecAXPY(e_P_F_ej, -ip(phi0[j], phi1[k]), phi1[k]); CHKERRQ(code); 

        code = VecAXPY(e_P_E_fi, -ip(phi1[i], phi0[k]), phi0[k]); CHKERRQ(code); 
        code = VecAXPY(e_P_E_fj, -ip(phi1[j], phi0[k]), phi0[k]); CHKERRQ(code); 
      }
      code = VecAXPY(e_P_F_ei, 1.0, phi0[i]);
      code = VecAXPY(e_P_F_ej, 1.0, phi0[j]);
      /* E */
      Aij = ip(e_P_F_ei, e_P_F_ej);
      code = MatSetValue(A, i, j, Aij, INSERT_VALUES); CHKERRQ(code);
      code = MatSetValue(A, j, i, Aij, INSERT_VALUES); CHKERRQ(code);
      /* F */
      Aij = ip(e_P_E_fi, e_P_E_fj);
      code = MatSetValue(A, n+i, n+j, Aij, INSERT_VALUES); CHKERRQ(code);
      code = MatSetValue(A, n+j, n+i, Aij, INSERT_VALUES); CHKERRQ(code);
      /* G*/
      Aij = ip(e_P_F_ei, e_P_E_fj);
      code = MatSetValue(A, n, n+j, Aij, INSERT_VALUES); CHKERRQ(code);
      code = MatSetValue(A, n+j, n, Aij, INSERT_VALUES); CHKERRQ(code);
    }
  }
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

  cout << "MAIN: solving generalized eigen-value problem." << endl;
  vector<vector<double>> e_values = la.eigenValueProblem(A, M, "GNHEP"); 
  cout << "MAIN: d(E,F) = " << max(e_values[0]) << endl; 

  ofstream file(par.dirResults() + "/olga" + to_string(par.nbModes()) + ".txt");
  file << max(e_values[0]) << endl;
  
 
  MatDestroy(&A);
  MatDestroy(&M);
  SlepcFinalize();
}
