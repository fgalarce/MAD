/*=============================================================================
 This file is part of MDMA.
 Multi-physics and Data assimilation for Medical Applications
 Copyright (C) 2020, 

  Felipe Galarce at INRIA.
 
 MDMA is free software; you can redistribute it and/or modify it under
 the terms of the GNU Lesser General Public License as published by the Free
 Software Foundation; version 2.1.

 MDMA is distributed in the hope that it will be useful, but WITHOUT ANY
 WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with MDMA. If not, see http://www.gnu.org/licenses/.
 =============================================================================*/

#include <ultra-4d-flow.hpp>

int main(int argc, char *argv[]){

  assert(argc == 2);

  string data_file = argv[1];
  
  /* Deploy Petsc */
  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  /* Initialize MDMA objects */
  Parameters par(data_file);
  par.print();

  IO io;
  io.initialize(par);

  Geometry geo;
  geo.initialize(par, io);

  CFD cfd;
  cfd.initialize(par, io, geo);

  InnerProduct ip;
  ip.initialize(par, geo);

  UsMeasures us;
  us.initialize(par, geo, ip);

  ofstream outfile(par.dirResults() + "/error.txt");

  int nbVertices = io.nbVertices();
  int nbDofs = 3 * io.nbVertices(); 
  int m = us.nbMeasures();
  int n = par.nbModes(); /* dim(Vn) */
  vector<vector<int>> N_mapping(5);
  N_mapping[0].push_back(80);
  N_mapping[0].push_back(75);
  N_mapping[0].push_back(70);
  N_mapping[0].push_back(65);
  N_mapping[0].push_back(90);

  N_mapping[1].push_back(85);
  N_mapping[1].push_back(90);
  N_mapping[1].push_back(80);
  N_mapping[1].push_back(75);
  N_mapping[1].push_back(100);

  N_mapping[2].push_back(65);
  N_mapping[2].push_back(100);
  N_mapping[2].push_back(90);
  N_mapping[2].push_back(105);
  N_mapping[2].push_back(95);

  N_mapping[3].push_back(60);
  N_mapping[3].push_back(75);
  N_mapping[3].push_back(75);
  N_mapping[3].push_back(80);
  N_mapping[3].push_back(75);

  N_mapping[4].push_back(55);
  N_mapping[4].push_back(55);
  N_mapping[4].push_back(60);
  N_mapping[4].push_back(75);
  N_mapping[4].push_back(65);

  ofstream mu_file(par.dirResults() + "/mu" + to_string(n) + ".txt");

  Mat W = mat(nbDofs + nbVertices, m, "dense");
  for (int i = 0; i < m; i++){
    vector<double> rr_vec(nbDofs + nbVertices);
    VecGetValues(us.rieszRepresenters()[i], nbDofs, &range(nbDofs)[0], &rr_vec[0]);
    MatSetValues(W, nbDofs, &range(nbDofs)[0], 1, &i, &rr_vec[0], INSERT_VALUES);
  }
  MatAssemblyBegin(W, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(W, MAT_FINAL_ASSEMBLY);
  Mat W_tr = transpose(W);

  /* Do computations for 25 windows */
  for (int idWinT = 0; idWinT < 5; idWinT++){
    for (int idWinHR = 0; idWinHR < 5; idWinHR++){
      cout << "--> window coordinates: " << idWinT << " " << idWinHR << endl;

      int N = N_mapping[idWinT][idWinHR];

      /* Load reduced model */
      vector<Vec> phi_u(N);
      vector<Vec> phi_p(N);
      vector<Vec> phi(N);
      for (int i = 0; i < N; i++){
        /* From ensight */
        vec(phi_u[i], nbDofs);
        vec(phi_p[i], nbVertices);
        loadVec(phi_u[i], par.dirModel() + "/window_time" + to_string(idWinT) + "_HR" + to_string(idWinHR) + "/mode_velocity." + wildcard(i) + ".vct");
        loadVec(phi_p[i], par.dirModel() + "/window_time" + to_string(idWinT) + "_HR" + to_string(idWinHR) +"/mode_pressure." + wildcard(i) + ".scl");
   
        vector<double> phi_u_vec(nbDofs); 
        vector<double> phi_p_vec(nbVertices); 
        code = VecGetValues(phi_u[i], nbDofs, &range(nbDofs)[0], &phi_u_vec[0]); CHKERR(code);
        code = VecGetValues(phi_p[i], nbVertices, &range(nbVertices)[0], &phi_p_vec[0]); CHKERR(code);
      
        /* To mixed petsc u-p */
        vec(phi[i], nbDofs + nbVertices);
        code = VecSetValues(phi[i], nbDofs, &range(nbDofs)[0], &phi_u_vec[0], INSERT_VALUES); CHKERR(code);
        code = VecSetValues(phi[i], nbVertices, &range(nbDofs, nbDofs + nbVertices)[0], &phi_p_vec[0], INSERT_VALUES); CHKERR(code);
        code = VecAssemblyBegin(phi[i]); CHKERR(code);
        code = VecAssemblyEnd(phi[i]); CHKERR(code);
      }

      cout << "--> assemble Tsi, an approximated basis for W^perp. Approximated dimension: " << N << endl;
      vector<Vec> Tsi(N);
      for (int i = 0; i < N; i++){
        cout << i << " ";
        vec(Tsi[i], nbDofs + nbVertices);
        VecCopy(phi[i], Tsi[i]);
        for (int j = 0; j < m; j++){
          code = VecAXPY(Tsi[i], -1.0*ip(us.rieszRepresenters()[j], phi[i]), us.rieszRepresenters()[j]); CHKERRQ(code);
        }
        for (int k = 0; k < i; k++){
          code = VecAXPY(Tsi[i], -1.0*ip(Tsi[k], Tsi[i]), Tsi[k]); CHKERRQ(code);
        }
      }
      cout << endl;

      cout << "--> Assemble P_Vn tsi_i" << endl;
      Mat M = mat(N,N, "dense");
      vector<Vec> proj_tsi_Vn(N);
      for (int i = 0; i < N; i++){
        cout << i << " ";
        zeros(proj_tsi_Vn[i], nbDofs + nbVertices);
        for (int j = 0; j < n; j++){
          code = VecAXPY(proj_tsi_Vn[i], ip(phi[j], Tsi[i]), phi[j]); CHKERRQ(code);
        }        
      }
      cout << endl;

      cout << "--> Assemble M = < tsi_i - P_Vn tsi_i , tsi_j - P_Vn tsi_j >" << endl;
      for (int i = 0; i < N; i++){
        cout << i << " ";
        double M_ij;
        Vec diff1 = vec(nbDofs + nbVertices);
        code = VecCopy(Tsi[i], diff1); CHKERR(code);
        code = VecAXPY(diff1, -1.0, proj_tsi_Vn[i]); CHKERRQ(code);
        for (int j = 0; j < N; j++){
          Vec diff2 = vec(nbDofs + nbVertices);      
          code = VecCopy(Tsi[j], diff2); CHKERRQ(code);
          code = VecAXPY(diff2, -1.0, proj_tsi_Vn[j]); CHKERRQ(code);
          code = MatSetValue(M, i, j, ip(diff1, diff2), INSERT_VALUES); CHKERRQ(code);
          VecDestroy(&diff2);
        }
        VecDestroy(&diff1);
      }
      code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
      code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
      cout << endl << "--> |M| = " << norm(M) << endl;;

      cout << "--> Assemble Q" << endl;
      Mat Q = mat(N,N, "dense");
      for (int i = 0; i < N; i++){
        cout << i << " ";
        double Q_i = cfd.integralSurface(Tsi[i], 3) / cfd.geo.computeBoundaryArea(3) - cfd.integralSurface(Tsi[i], 2) / cfd.geo.computeBoundaryArea(2);
        for (int j = 0; j < N; j++){
          double Q_ij = Q_i * cfd.integralSurface(Tsi[j], 3) / cfd.geo.computeBoundaryArea(3) - cfd.integralSurface(Tsi[j], 2) / cfd.geo.computeBoundaryArea(2);
          code = MatSetValue(Q, i, j, Q_ij, INSERT_VALUES); CHKERRQ(code);
        }
      }
      code = MatAssemblyBegin(Q, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
      code = MatAssemblyEnd(Q, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
      cout << endl << "--> |Q| = " << norm(Q) << endl;;

      cout << "--> Solving generalized eigen-value problem Q c = a M c" << endl;
      EPS eps;
      code = EPSCreate(PETSC_COMM_WORLD, &eps); CHKERRQ(code);
      code = EPSSetOperators(eps, Q, M); CHKERRQ(code);
      code = EPSSetProblemType(eps, EPS_GNHEP); CHKERRQ(code);
      code = EPSSetType(eps, EPSKRYLOVSCHUR); CHKERRQ(code);
      code = EPSSetDimensions(eps, N, PETSC_DECIDE, PETSC_DECIDE); CHKERRQ(code);
      code = EPSSetFromOptions(eps); CHKERRQ(code);
      code = EPSSolve(eps); CHKERRQ(code);
      double maxev = 0.0;
      for (int i = 0; i < N; i++){
        double eigr, eigi;
        code = EPSGetEigenpair(eps, i, &eigr, &eigi, NULL, NULL);
        if (sqrt(eigr*eigr + eigi*eigi) > maxev){
          maxev = sqrt(eigr*eigr + eigi*eigi);
        }
      }
      cout << "--> mu = " << sqrt(maxev) << endl;
      mu_file << sqrt(maxev) << endl;
      EPSDestroy(&eps);
      MatDestroy(&Q);
      MatDestroy(&M);
      for (int i = 0; i < N; i++){
        VecDestroy(&phi_u[i]);
        VecDestroy(&phi_p[i]);
        VecDestroy(&phi[i]);
        VecDestroy(&Tsi[i]);
        VecDestroy(&proj_tsi_Vn[i]);
      }
    }
  }

  mu_file.close();
  SlepcFinalize();

}
