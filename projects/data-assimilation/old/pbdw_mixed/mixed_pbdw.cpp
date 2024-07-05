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

  Model Vn;
  Vn.initialize(par, io.nbVertices());

  Geometry geo;
  geo.initialize(par, io);

  InnerProduct ip;
  ip.initialize(par, geo);

  UsMeasures us;
  us.initialize(par, geo, ip);

  int nbVertices = io.nbVertices();
  int nbDofs = 4 * io.nbVertices(); 
  int nbDofsVelocity = 3 * io.nbVertices(); 
  int nbDofsPressure = io.nbVertices(); 
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
  error_least_squares << "|up - up*|_H1xL2 ";
  error_least_squares << "|up|_H1xL2 ";
  error_least_squares << "|up - P_Vn up|_H1xL2 ";
  error_least_squares << "|u - u*|_H1 ";
  error_least_squares << "|u|_H1 ";
  error_least_squares << "|p - p*|_L2 ";
  error_least_squares << "|p|_L2 ";
  error_least_squares << "|grad(u) - grad(u*)|_L2 ";
  error_least_squares << "|grad(u)|_L2" << endl;
  ofstream error_pbdw(par.dirResults() + "/error_pbdw.txt");
  error_pbdw << "|up - up*|_H1xL2 ";
  error_pbdw << "|up|_H1xL2 ";
  error_pbdw << "|up - P_Vn up|_H1xL2 ";
  error_pbdw << "|u - u*|_H1 ";
  error_pbdw << "|u|_H1 ";
  error_pbdw << "|p - p*|_L2 ";
  error_pbdw << "|p|_L2 ";
  error_pbdw << "|grad(u) - grad(u*)|_L2 ";
  error_pbdw << "|grad(u)|_L2" << endl;
  ofstream modalContributions(par.dirResults() + "/modal_contributions.txt");

  for (double time : range(par.UStimeStart(), par.UStimeEnd(), par.UStimeStep())){

    int iteration;
    if ( abs(time - floor(time)) > 0.5){
      iteration = floor(time/0.004) + 1;
    } else {
      iteration = floor(time/0.004);
    }

    cout << "\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =" << endl;
    cout << "Reconstruction time: " << time << endl;
    cout << "Normalized time: " << time*par.heartRate()/60  << endl;
    cout << "Iteration: " << iteration << endl;

    /* Compute Doppler from perfect model */
    us.computeSynthetic(time);

    /* Get basis of Vn and it's dimension for current time window */
    cout << "Curent time window: " << Vn.indexCurrentBasis() << "/" << par.nbTimeWindows() << endl;
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

          loadMat(G, "../../data/model/" + par.method() + "_linear_H1/window_time" + to_string(Vn.indexCurrentBasis()) 
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
          vector<int> index_rows = range(m);
          for (PetscInt j = 0; j < n; j++){
            cout << "Assembling G: " << float(j)/float(n)*100 << " %      \r";
            for (size_t i = 0; i < m; i++){
              Gij[i] = ip(phi[j], us.rieszRepresenters()[i]);
            }
            code = MatSetValues(G, m, &index_rows[0], 1, &j, &Gij[0], INSERT_VALUES); CHKERRQ(code); 
            code = MatSetValues(G_tr, 1, &j, m, &index_rows[0], &Gij[0], INSERT_VALUES); CHKERRQ(code); 
          }

          code = MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
          code = MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
          code = MatAssemblyBegin(G_tr, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
          code = MatAssemblyEnd(G_tr, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 

          vector<double> sValuesG = getSingularValues(G, n);
          exportData(sValuesG, par.dirResults() + "/G_sv_" + wildcard( Vn.indexCurrentBasis()) + ".txt");
          cout << "Singular values of G: " << sValuesG << endl;
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

    Vec velocity = vec(nbDofsVelocity);
    Vec pressure = vec(nbDofsPressure) ;

    /* Export reconstruction to Ensight */
    double velocity_array[nbVertices*3];
    double pressure_array[nbVertices];

    code = VecGetValues(v_star, nbDofsVelocity, &range(nbDofsVelocity)[0], velocity_array); CHKERRQ(code); 
    code = VecGetValues(v_star, nbDofsPressure, &range(nbDofsVelocity, nbDofsVelocity + nbDofsPressure)[0], pressure_array); CHKERRQ(code); 

    code = VecSetValues(velocity, 3*nbVertices, &range(nbDofsVelocity)[0], velocity_array, INSERT_VALUES); CHKERRQ(code);
    code = VecAssemblyBegin(velocity); CHKERRQ(code);
    code = VecAssemblyEnd(velocity); CHKERRQ(code);

    code = VecSetValues(pressure, nbVertices, &range(nbDofsPressure)[0], pressure_array, INSERT_VALUES); CHKERRQ(code);
    code = VecAssemblyBegin(pressure); CHKERRQ(code);
    code = VecAssemblyEnd(pressure); CHKERRQ(code);

    if (!par.minimalOutput()) io.writeState(velocity, "v_star", time);
    if (!par.minimalOutput()) io.writeState(pressure, "p_star", time); 

    /*
        Load Ground Truth
    */
    Vec u = vec(nbDofsVelocity);
    loadVec(u, par.dirSyntheticField() + "velocity." + wildcard(iteration) +  ".vct");

    Vec target;
    PetscReal u_array[nbVertices*3];
    PetscReal p_array[nbVertices];

    Vec p = vec(io.nbVertices());
    code = VecCreate(PETSC_COMM_WORLD, &target); CHKERRQ(code);

    code = VecSetSizes(target, PETSC_DECIDE, nbVertices*4); CHKERRQ(code);
    code = VecSetFromOptions(target); CHKERRQ(code);

    loadVec(p, par.dirSyntheticField() + "pressure." + wildcard(iteration) + ".scl");

    code = VecGetValues(u, 3*nbVertices, &range(nbDofsVelocity)[0], u_array); CHKERRQ(code); 
    code = VecGetValues(p, nbVertices, &range(nbDofsPressure)[0], p_array); CHKERRQ(code); 

    code = VecSetValues(target, 3*nbVertices, &range(nbDofsVelocity)[0], u_array, INSERT_VALUES); CHKERRQ(code);
    code = VecSetValues(target, nbVertices, &range(nbDofsVelocity, nbDofsVelocity + nbDofsPressure)[0], p_array, INSERT_VALUES); CHKERRQ(code);

    code = VecAssemblyBegin(target); CHKERRQ(code);
    code = VecAssemblyEnd(target); CHKERRQ(code);


    if (!par.minimalOutput()) io.writeState(u, "velocity", time);

    /*
        Compute error least-squares 
    */
   
    /* Compute model error */
    Vec model_error = vec(nbDofs);
    VecCopy(target, model_error);
    for (int i = 0; i < n; i++){
      double proj_Uvn;
      proj_Uvn = ip(phi[i], target); 
      code = VecAXPY(model_error,  -proj_Uvn, phi[i]); CHKERRQ(code);
    }
    double normModelError = sqrt(ip(model_error, model_error));

    /* Compute relative error */
    Vec error_x = zeros(nbDofs); 
    Vec error_u = zeros(3*nbVertices); 
    Vec error_p = zeros(nbVertices); 
    double normU, normP, normX, norm_gradU; 
    double errorP, errorU, error_gradU;

    code = VecAXPY(error_x,  1.0, target); CHKERRQ(code);
    code = VecAXPY(error_x, -1.0, v_star); CHKERRQ(code);
    if (!par.minimalOutput()) io.writeState(error_x, "error", time);

    normX = sqrt(ip(target, target));
    double errorX = sqrt(ip(error_x, error_x));

    code = VecAXPY(error_u,  1.0, velocity); CHKERRQ(code);
    code = VecAXPY(error_u, -1.0, u); CHKERRQ(code);
    code = VecAXPY(error_p,  1.0, pressure); CHKERRQ(code);
    code = VecAXPY(error_p, -1.0, p); CHKERRQ(code);
    errorP = sqrt(ip(error_p, error_p));
    errorU = sqrt(ip(error_u, error_u));
    error_gradU = sqrt(ip.only_grad(error_u, error_u));
    normU = sqrt(ip(u, u));
    normP = sqrt(ip(p, p));
    norm_gradU = sqrt(ip.only_grad(u, u));

    cout << "  Summary for least-squares up:" << endl;
    cout << scientific << "  | up | = " << norm(target) << endl;
    cout << scientific << "  | up* | = " << norm(v_star) << endl;
    cout << scientific << "  | up - P_vn up | / | up | = " <<  normModelError/normX<< endl;
    cout << scientific << "  | up - up* | / |up| = " << errorX/normX << endl;
    cout << scientific << "  | u - u* | / | u | = " << errorU/normU << endl;
    cout << scientific << "  | grad_u - grad_u* | / | grad_u | = " << error_gradU/norm_gradU << endl;
    cout << scientific << "  | p - p* | / | p | = " << errorP/normP << endl;

    error_least_squares << scientific << errorX << " ";
    error_least_squares << scientific << normX << " ";
    error_least_squares << scientific << normModelError << " ";
    error_least_squares << scientific << errorU << " ";
    error_least_squares << scientific << normU << " ";
    error_least_squares << scientific << errorP << " ";
    error_least_squares << scientific << normP << " ";
    error_least_squares << scientific << error_gradU << " ";
    error_least_squares << scientific << norm_gradU << endl;

    /*
      Compute error pbdw 
    */

    /* residual */
    Vec res = vec(m);
    Vec Gc = vec(m);
    MatMult(G, x, Gc);
    VecCopy(us.measures(), res);
    VecAXPY(res, -1.0, Gc);

    for (PetscInt j = 0; j < m; j++){
      double rj;
      VecGetValues(res, 1, &j, &rj);
      VecAXPY(v_star, rj, us.rieszRepresenters()[j]);
    }

    /* Compute relative error */
    error_x = zeros(nbDofs); 
    error_u = zeros(3*nbVertices); 
    error_p = zeros(nbVertices); 

    code = VecAXPY(error_x,  1.0, target); CHKERRQ(code);
    code = VecAXPY(error_x, -1.0, v_star); CHKERRQ(code);
    if (!par.minimalOutput()) io.writeState(error_x, "error", time);

    code = VecAXPY(error_u,  1.0, velocity); CHKERRQ(code);
    code = VecAXPY(error_u, -1.0, u); CHKERRQ(code);
    code = VecAXPY(error_p,  1.0, pressure); CHKERRQ(code);
    code = VecAXPY(error_p, -1.0, p); CHKERRQ(code);

    errorX = sqrt(ip(error_x, error_x));
    errorP = sqrt(ip(error_p, error_p));
    errorU = sqrt(ip(error_u, error_u));
    error_gradU = sqrt(ip.only_grad(error_u, error_u));

    cout << "  Summary for pbdw up:" << endl;
    cout << scientific << "  | up | = " << norm(target) << endl;
    cout << scientific << "  | up* | = " << norm(v_star) << endl;
    cout << scientific << "  | up - P_vn up | / | up | = " <<  normModelError/normX<< endl;
    cout << scientific << "  | up - up* | / |up| = " << errorX/normX << endl;
    cout << scientific << "  | u - u* | / | u | = " << errorU/normU << endl;
    cout << scientific << "  | grad_u - grad_u* | / | grad_u | = " << error_gradU/norm_gradU << endl;
    cout << scientific << "  | p - p* | / | p | = " << errorP/normP << endl;

    error_pbdw << scientific << errorX << " ";
    error_pbdw << scientific << normX << " ";
    error_pbdw << scientific << normModelError << " ";
    error_pbdw << scientific << errorU << " ";
    error_pbdw << scientific << normU << " ";
    error_pbdw << scientific << errorP << " ";
    error_pbdw << scientific << normP << " ";
    error_pbdw << scientific << error_gradU << " ";
    error_pbdw << scientific << norm_gradU << endl;

    /* update time */
    time = time + par.UStimeStep();
   
     
    code = VecDestroy(&error_x); CHKERRQ(code);
    code = VecDestroy(&error_u); CHKERRQ(code);
    code = VecDestroy(&error_p); CHKERRQ(code);
    code = VecDestroy(&v_star); CHKERRQ(code);
    code = VecDestroy(&u); CHKERRQ(code);
    code = VecDestroy(&model_error); CHKERRQ(code);
    code = VecDestroy(&target); CHKERRQ(code);
    code = VecDestroy(&velocity); CHKERRQ(code);
    code = VecDestroy(&pressure); CHKERRQ(code);
    
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
