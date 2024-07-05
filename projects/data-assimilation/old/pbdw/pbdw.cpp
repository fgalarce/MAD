#include <ultra-4d-flow.hpp>


int main(int argc, char *argv[]){

  assert(argc == 2);

  string data_file = argv[1];
  
  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

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

  UsMeasures us;
  us.initialize(par, geo, la);

  Model Vn;
  Vn.initialize(par, geo, la);

  int nbVertices = io.nbVertices();
  int nbDofs = 3 * io.nbVertices(); 
  int m = us.nbMeasures();
  int n; /* dim(Vn) */

  /* Declare Cross-Gramian matrix and normal equation A x = b objects */
  Mat G, G_tr, A;
  Vec b, x;

  if (!par.adaptive_nbModes()){
    n = par.nbModes();
    mat(G, m, n);
    mat(G_tr, n, m);
    mat(A, n, n);
    vec(b, n);
    if (!par.restrainedLeastSquares())
      vec(x,n);
  }

  /* Files to write error and coefficients in reduced model */
  ofstream error_least_squares(par.dirResults() + "/error_ls.txt");
  error_least_squares << "error_u\tnormU\terror_model\terror_gradU\tnorm_gradU" << endl;
  ofstream error_pbdw(par.dirResults() + "/error_pbdw.txt");
  error_pbdw << "error_u\tnormU\terror_model\terror_gradU\tnorm_gradU" << endl;
  ofstream modalContributions(par.dirResults() + "/modal_contributions.txt");

  for (double time : range(par.UStimeStart(), par.UStimeEnd(), par.UStimeStep())){

    int iteration;
    if ( abs(time - floor(time)) > 0.5){
      iteration = floor(time/par.CFDtimeStep()) + 1;
    } else {
      iteration = floor(time/par.CFDtimeStep());
    }

    cout << "\n- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
    cout << "Reconstruction time: " << time << endl;
    cout << "Normalized time: " << time*par.heartRate()/60  << endl;
    cout << "Iteration: " << iteration << endl;

    /* Compute Doppler from perfect model */
    us.computeSynthetic(iteration);

    /* Get basis of Vn and it's dimension for current time window */
    cout << "Curent time window: " << Vn.indexCurrentBasis() << "/" << par.nbTimeWindows() - 1 << endl;
    vector<Vec> phi = Vn.basis(time);
    n = Vn.nbModesMapping()[Vn.indexCurrentBasis()];

    if (Vn.updateBasis()){

      /* Either assemble or load Cross-Gramian matrix*/
      if (par.adaptive_nbModes()){
        cout << "Current number of modes: " << n << endl;
        vec(b,    n);
        vec(x,    n);
        mat(G,    m, n);
        mat(G_tr, n, m);
        mat(A,    n, n);
        loadMat(G, par.dirModel() + "/" + par.method() + "_" + par.type()  + "_"+ par.innerProduct() 
                  + "/window_time" + to_string(Vn.indexCurrentBasis()) + "_HR" + to_string( Vn.coord_HR())  
                  + "/G_nstar_" + to_string(n) + ".bin");    

        for (PetscInt j = 0; j < n; j++){
          cout << "Assembling G_tr: " << float(j)/float(n)*100 << " %      \r";
          for (PetscInt i = 0; i < m; i++){
            vector<double> Gij_(1);
            MatGetValues(G, 1, &i, 1, &j, &Gij_[0]);
            MatSetValue(G_tr, j, i, Gij_[0], INSERT_VALUES);
          }
        }
        code = MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
        code = MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
        code = MatAssemblyBegin(G_tr, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
        code = MatAssemblyEnd(G_tr, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 

      } else {

        if (par.loadMeasuresAndG()){

          loadMat(G, "../../data/model/pod_linear_H1/window_time" + to_string(Vn.indexCurrentBasis())
                     + "_HR" + to_string(Vn.coord_HR()) + "/G_" + to_string(n) + "modes.bin");

          cout << "Assembling G_tr row: [ ";
          for (PetscInt i = 0; i < m; i++){
            cout << i << " ";
            for (PetscInt j = 0; j < n; j++){
              vector<double> Gij_(1);
              code = MatGetValues(G, 1, &i, 1, &j, &Gij_[0]); CHKERRQ(code);
              code = MatSetValue(G_tr, j, i, Gij_[0], INSERT_VALUES); CHKERRQ(code);
            }
          }
          cout << " ]" << endl;
          code = MatAssemblyBegin(G_tr, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
          code = MatAssemblyEnd(G_tr, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 

        } else {

          vector<double> Gij(m);
          for (PetscInt j = 0; j < n; j++){
            cout << "Assembling G: row " << j << "." << endl;
            for (size_t i = 0; i < m; i++){
              Gij[i] = ip(phi[j], us.rieszRepresenters()[i]);
            }
            code = MatSetValues(G, m, &range(m)[0], 1, &j, &Gij[0], INSERT_VALUES); CHKERRQ(code); 
            code = MatSetValues(G_tr, 1, &j, m, &range(m)[0], &Gij[0], INSERT_VALUES); CHKERRQ(code); 
          }

          code = MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
          code = MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
          code = MatAssemblyBegin(G_tr, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
          code = MatAssemblyEnd(G_tr, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 

          vector<double> sValuesG = getSingularValues(G, n);
          exportData(sValuesG, par.dirResults() + "/G_sv_" + wildcard( Vn.indexCurrentBasis()) + ".txt");
          cout << "PBDW: Singular values of G: " << endl << sValuesG << endl;
        }
      }
      code = MatMatMult(G_tr, G, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &A); CHKERRQ(code);
    }

    /* G^T l */
    code = MatMult(G_tr, us.measures(), b); CHKERRQ(code);

    if (par.restrainedLeastSquares()){ 
      /* - - - - - - - - - - - -
          Solve KKT conditions for constrained problem
         - - - - - - - - - - - - */
      /* 
          [ A = G^T G  | I | -I ] [ c   ] = [G^T y] 
          [  I         | 0 |  0 ] [\mu_1] = [\delta]
          [ -I         | 0 |  0 ] [\mu_2] = [\delta]
      */

      cout << "Solving normal equations without constraints." << endl;
      KSP ksp;
      code = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(code); 
      code = KSPSetOperators(ksp, A, A); CHKERRQ(code); 
      code = KSPSetType(ksp, KSPGMRES); CHKERRQ(code); 
      PC pc;
      code = KSPGetPC(ksp, &pc); CHKERRQ(code); 
      code = PCSetType(pc, PCASM); CHKERRQ(code); 
      code = KSPSetFromOptions(ksp); CHKERRQ(code); 
      vec(x, n);
      code = KSPSolve(ksp, b, x); CHKERRQ(code); 
      code = KSPDestroy(&ksp); CHKERRQ(code); 

      cout << "**********************" << endl;
      cout << "Evaluating duality constraints. " << endl;
      /* assemble delta */
      vector<double> delta_stl = Vn.delta(); /* set initial delta */
      vector<double> x_stl(n);
      VecGetValues(x, n, &range(n)[0], &x_stl[0]);
      for (int i : range(n)){
        if (x_stl[i] > delta_stl[i]){
          cout << "\tCondition " << i << " is active." << endl;
          x_stl[i] = delta_stl[i];
        } else if (x_stl[i] < -delta_stl[i]){
          cout << "\tCondition " << i << " is active." << endl;
          x_stl[i] = -delta_stl[i];
        }
      }
      VecSetValues(x, n, &range(n)[0], &x_stl[0], INSERT_VALUES);
      VecAssemblyBegin(x);
      VecAssemblyEnd(x);

    } else {
      cout << "Solving normal equations." << endl;
      KSP ksp;
      code = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(code); 
      code = KSPSetOperators(ksp, A, A); CHKERRQ(code); 
      code = KSPSetType(ksp, KSPGMRES); CHKERRQ(code); 
      PC pc;
      code = KSPGetPC(ksp, &pc); CHKERRQ(code); 
      code = PCSetType(pc, PCASM); CHKERRQ(code); 
      code = KSPSetFromOptions(ksp); CHKERRQ(code); 
      code = KSPSolve(ksp, b, x); CHKERRQ(code); 
      code = KSPDestroy(&ksp); CHKERRQ(code); 
    }

    /* Computing linear combination of functions in Vn */
    Vec v_star = zeros(nbDofs);
    double modal_contribution[n];
    code = VecGetValues(x, n, &range(n)[0], modal_contribution); CHKERRQ(code);
    cout << "Modal contributions: " << endl;
    for (int i = 0; i < n; i++){
      modalContributions << scientific << modal_contribution << " ";
      VecAXPY(v_star, modal_contribution[i], phi[i]);
      cout << modal_contribution[i] << " ";
    }
    modalContributions << endl;
    cout << endl;
    
//    VecView(b, PETSC_VIEWER_STDOUT_WORLD);
//    MatView(A, PETSC_VIEWER_STDOUT_SELF);
//    VecView(x, PETSC_VIEWER_STDOUT_SELF);
//    exit(1);

    /* Export reconstruction to Ensight */
    if (!par.minimalOutput()) io.writeState(v_star, "v_star", time);

    /* Compute error in v_star */
    Vec u = vec(nbDofs);
    loadVec(u, par.dirSyntheticField() + "phi_vector." + wildcard(iteration) +  par.modelFormat());
    if (!par.minimalOutput()) io.writeState(u, "velocity", time);
    double normU = sqrt(ip(u,u));

    Vec error = vec(nbDofs);
    VecCopy(u, error);
    VecAXPY(error, -1.0, v_star);
    if (!par.minimalOutput()) io.writeState(error, "error_ls", time);

    double normErrorU = sqrt(ip(error,error));
    double normGradU;
    double normErrorGradU;
    if (par.innerProduct() == "H1") normGradU = sqrt(ip.only_grad(u,u)); 
    if (par.innerProduct() == "H1") normErrorGradU = sqrt(ip.only_grad(error,error));

    Vec model_error = vec(nbDofs);
    VecCopy(u, model_error);
    for (int i = 0; i < n; i++){
      PetscScalar proj_Uvn;
      proj_Uvn = ip(phi[i], u); 
      code = VecAXPY(model_error,  -proj_Uvn, phi[i]); CHKERRQ(code);
    }
    double normErrorModel = sqrt(ip(model_error, model_error));

    if (par.innerProduct() == "H1") error_least_squares << scientific << normErrorU << " " << normU << " " << normErrorModel << " " << normErrorGradU << " " << normGradU << endl;
    if (par.innerProduct() == "H1") error_least_squares << scientific << normErrorU << " " << normU << " " << normErrorModel << endl;
    
    cout << "  Summary for least-squares: " << endl;    
    cout << scientific << "\t| u | = " << normU << endl;
    cout << scientific << "\t| u* | = " << sqrt(ip(v_star,v_star)) << endl;
    cout << scientific << "\t| u - P_vn u | / | u | = " <<  normErrorModel / normU<< endl;
    cout << scientific << "\t| u - u* | / |u| = " << normErrorU / normU << endl;
    if (par.innerProduct() == "H1") cout << scientific << "\t| grad_u - grad_u* | / | grad_u | = " << normErrorGradU / normGradU << endl;

    /*
        Compute PBDW field 
    */

    /* compute residual r = l - Gc*/
    Vec res = vec(m);
    Vec Gc = vec(m);
    code = MatMult(G, x, Gc); CHKERRQ(code);
    VecCopy(us.measures(), res);
    VecAXPY(res, -1.0, Gc);

    for (PetscInt j = 0; j < m; j++){
      double rj;
      VecGetValues(res, 1, &j, &rj);
      VecAXPY(v_star, rj, us.rieszRepresenters()[j]);
    }

    /* Compute error in pbdw */
    error = zeros(nbDofs);
    VecCopy(u, error);
    VecAXPY(error, -1.0, v_star);
    if (!par.minimalOutput()) io.writeState(error, "error_pbdw", time);

    normErrorU = sqrt(ip(error,error));
    if (par.innerProduct() == "H1") normGradU = sqrt(ip.only_grad(u,u)); 
    if (par.innerProduct() == "H1") normErrorGradU = sqrt(ip.only_grad(error,error));

    model_error = zeros(nbDofs);
    VecCopy(u, model_error);
    for (int i = 0; i < n; i++){
      PetscScalar proj_Uvn;
      proj_Uvn = ip(phi[i], u); 
      code = VecAXPY(model_error,  -proj_Uvn, phi[i]); CHKERRQ(code);
    }
    normErrorModel = sqrt(ip(model_error, model_error));

    error_pbdw << scientific << normErrorU << " " << normU << " " << normErrorModel << " " << normErrorGradU << " " << normGradU << endl;
    cout << "  Summary for PBDW: " << endl;    
    cout << scientific << "\t| u | = " << normU << endl;
    cout << scientific << "\t| u* | = " << sqrt(ip(v_star,v_star)) << endl;
    cout << scientific << "\t| u - P_vn u | / | u | = " <<  normErrorModel / normU<< endl;
    cout << scientific << "\t| u - u* | / |u| = " << normErrorU / normU << endl;
    if (par.innerProduct() == "H1") cout << scientific << "\t| grad_u - grad_u* | / | grad_u | = " << normErrorGradU / normGradU << endl;

    /* update time */
    time = time + par.UStimeStep();
    
    code = VecDestroy(&error); CHKERRQ(code);
    code = VecDestroy(&model_error); CHKERRQ(code);
    code = VecDestroy(&v_star); CHKERRQ(code);
    code = VecDestroy(&u); CHKERRQ(code);
  }
  
  /* Free the ram memory */
  code = MatDestroy(&G); CHKERRQ(code); 
  code = MatDestroy(&A); CHKERRQ(code); 
  code = MatDestroy(&G_tr); CHKERRQ(code); 
  if (!par.adaptive_nbModes()){
    code = VecDestroy(&b); CHKERRQ(code);
    code = VecDestroy(&x); CHKERRQ(code);
  }

  modalContributions.close();
  error_least_squares.close();
  error_pbdw.close();
  us.finalize();
  Vn.finalize();
  SlepcFinalize();

}
