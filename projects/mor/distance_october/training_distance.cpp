#include<ultra-4d-flow.hpp>

#define EPSI 0.00001

int main(int argc, char *argv[]){

  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  /* Parse data file */
  assert(argc == 4);
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
  vector<Vec> phi(n);
  vector<Vec> phi_template(n);
  for(int i=0;i<n;i++){
    vec(phi[i], nbDofs);
    vec(phi_template[i], nbDofs);
    if (stoi(argv[2]) == stoi(argv[3])){
      loadVec(phi[i], par.model() + "/mode." + wildcard(i) + ".vct");
    } else {
      loadVec(phi[i], par.model() + "/I_Tpod." + wildcard(i) + ".vct");
    }
    loadVec(phi_template[i], par.templateModel() + "/mode." + wildcard(i) + ".vct");
  }

  if (!la.checkOrthonormality(phi)){
    cout << "MAIN: basis from " << par.model() << " are not orthonormal." << endl;
    exit(1);
  } 

  if (!la.checkOrthonormality(phi_template)){
    cout << "MAIN: basis from " << par.templateModel() << " are not orthonormal." << endl;
    exit(1);
  } 

  double eij = 0.0;
  int snaps_per_sim = 25;
  for (int k = 0; k < par.nbSimulations(); k++){
    for (int snapId = 0; snapId < snaps_per_sim; snapId++){
      Vec e = zeros(nbDofs);
      Vec u = vec(nbDofs);
      Vec Pu = zeros(nbDofs);
      loadVec(u, par.maniFolder() + "sim" + wildcard(k) + "/phi_vector." + wildcard(snapId) + ".vct");
      for(int m=0; m<n; m++){
        VecAXPY(e, ip(u, phi[m]), phi[m]);
        VecAXPY(Pu, ip(u, phi_template[m]), phi_template[m]);
      }
      VecAXPY(e, -1.0, Pu);
      if (sqrt(ip(u,u)) > EPSI) {
        eij = eij + sqrt(ip(e,e) / ip(Pu,Pu));
      }
    }
  }
  eij = eij / ((double)par.nbSimulations() + (double)snaps_per_sim);

  cout << "MAIN: e_ij = " << eij << endl;
  ofstream eij_file(par.dirResults() + "/e" + argv[2] + "_" + argv[3]);
  eij_file << eij;
  eij_file.close();
  SlepcFinalize();
  
  return 0;
}
