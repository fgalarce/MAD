#include<ultra-4d-flow.hpp>

int main(int argc, char *argv[]){

  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  /* Parse data file */
  assert(argc == 3);
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
    loadVec(phi0[i], par.templateModel() + "/mode." + wildcard(i)+".vct");
    vec(phi1[i], nbDofs);
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

  Mat Phi0 = buildProjector(phi0);
  Mat Phi1 = buildProjector(phi1);
  Mat Phi0_tr = mat(n, nbDofs, "dense");
  code = MatTranspose(Phi0, MAT_INITIAL_MATRIX, &Phi0_tr); CHKERRQ(code);
  Mat Phi0_tr_M = mat(n, nbDofs);
  code = MatMatMult(Phi0_tr, ip.matrix, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Phi0_tr_M); CHKERRQ(code);
  Mat A = mat(n,n,"dense");
  code = MatMatMult(Phi0_tr_M, Phi1, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &A); CHKERRQ(code);

  cout << "MAIN: Norm E: " << norm(Phi0) << endl;
  cout << "MAIN: Norm F: " << norm(Phi1) << endl;
  cout << "MAIN: Norm F^T E: " << norm(Phi1) << endl;

  cout << "MAIN: Computing SVD(F^T (M) E)." << endl;
//  int l = n;
//  if (n > N){
//    l = N;
//  }
  vector<double> svalues = getSingularValues(A, n);
  cout << "MAIN: Computed " << svalues.size() << " singular values." << endl;
  vector<double> theta = arccos(svalues); 

  /* Asimov */
  double dA = theta[par.nbModes()-1];

  /* Binet-Cauchy */
  double dBC = 1;
  for (int i=0; i<par.nbModes(); i++){
    dBC = dBC*cos(theta[i])*cos(theta[i]);
    cout << dBC << endl;
  }
  dBC = sqrt(1.0 - dBC);

  /* Chordal */
  double dC = 0;
  for (int i=0; i<par.nbModes(); i++){
    dC += sin(theta[i])*sin(theta[i]);
  }
  dC = sqrt(dC);

  /* Fubini-Study */
  double dFS = 1;
  for (int i=0; i<par.nbModes(); i++){
    dFS = dFS*svalues[i];
  }
  dFS = acos(dFS);

  /* Martin */
  double dM = 1;
  for (int i=0; i<par.nbModes(); i++){
    dM = dM/(svalues[i]*svalues[i]);
  }
  dM = sqrt(log(dM));

  /* Procrustes */
  double dP = 1;
  for (int i=0; i<par.nbModes(); i++){
    dP += dP * sin(theta[i]/2.0)*sin(theta[i]/2.0);
  }
  dP = 2*sqrt(dP);

  double dProj = sin(theta[par.nbModes()-1]);
  double dSpec = 2*sin(theta[par.nbModes()-1]/2);

  double dGras = 0.0;
  for (int i=0; i<par.nbModes(); i++){
    dGras += theta[i] * theta[i];
  } 
  dGras = sqrt(dGras);
  
  double dH = sqrt(1 - min(svalues)*min(svalues));

  cout << "MAIN: S(A) = \n" << svalues << endl; 
  cout << "Main: theta(A) = \n" << theta << endl; 
  cout << endl; 
  cout << "MAIN: - - - - - SUMMARY - - - - -" << endl;
  cout << "MAIN: Asimov = " << dA   << endl;
  cout << "MAIN: Binet-Cauchy = " << dBC << endl;
  cout << "MAIN: Chordal = " << dC << endl;
  cout << "MAIN: Fubini-Study = " << dFS << endl;
  cout << "MAIN: Martin = " << dM << endl;
  cout << "MAIN: Procrustes = " << dP << endl;
  cout << "MAIN: Projection = " << dProj << endl;
  cout << "MAIN: Spectral = " << dSpec << endl;
  cout << "MAIN: Grassmanian = " << dGras << endl;
  cout << "MAIN: Hausdorff = " << dH << endl;

  ofstream file(par.dirResults() + "/distBC" + argv[2] + ".txt");
  file << dA << endl;
  file << dBC << endl;
  file << dC << endl;
  file << dFS << endl;
  file << dM << endl;
  file << dP << endl;
  file << dProj << endl;
  file << dSpec << endl;
  file << dGras << endl;
  file << dH << endl;
  file.close();

  MatDestroy(&Phi1);
  MatDestroy(&Phi0);
  MatDestroy(&Phi0_tr);
  MatDestroy(&A);
  SlepcFinalize();
}
