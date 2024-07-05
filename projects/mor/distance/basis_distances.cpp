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
    loadVec(phi0[i], par.model() + "/mode." + wildcard(i)+".vct");
    loadVec(phi1[i], par.templateModel() + "/I_Tpod." + wildcard(i) + ".vct");
  }

  if (!la.checkOrthonormality(phi1)){
    cout << "MAIN: basis from " << par.templateModel() << " are not orthonormal." << endl;
    exit(1);
  } 
  if (!la.checkOrthonormality(phi0)){
    cout << "MAIN: basis from " << par.model() << " are not orthonormal." << endl;
    exit(1);
  }

  for (int i = 0; i < par.nbModes(); i++){
    Vec error = vec(nbDofs);
    VecCopy(phi0[i], error);
    VecAXPY(error, -1.0, phi1[i]);
//    cout << "MAIN: distance = " << sqrt(ip(error,error)/ip(phi0[i], phi0[i])) << endl;
    cout << "MAIN: distance = " << sqrt(ip(error,error)) << endl;
  }

  SlepcFinalize();
}
