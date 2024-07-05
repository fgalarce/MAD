
#include<ultra-4d-flow.hpp>

int main(int argc, char *argv[]){

  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  /* Parse data file */
  assert(argc == 4);
  Parameters par(argv[1]);

  /* Initialize MDMA objects */
  par.print();


  IO io;
  io.initialize(par);

  Geometry geo;
  geo.initialize(par, io);

  int nbDofs = 3*geo.nbVertices;

  /* Compute work done by LDDMM field */
  CFD cfd;
  cfd.initialize(par, geo); 

  /* Load surface registration mappings */
  Vec lddmm_kl = vec(nbDofs);
  loadVec(lddmm_kl, par.surfaceMapping());

  Vec lddmm_lk = vec(nbDofs);
  loadVec(lddmm_lk, par.surfaceMappingInverse());

  /* T^T Mb T */
  cfd.assembleMassBoundary();
  Vec mtv = vec(nbDofs);
  code = MatMult(cfd.massBD(), lddmm_kl, mtv); CHKERRQ(code);
  double dij1;
  code = VecDot(mtv, lddmm_kl, &dij1); CHKERRQ(code);

  mtv = zeros(nbDofs);
  code = MatMult(cfd.massBD(), lddmm_lk, mtv); CHKERRQ(code);
  double dij2;
  code = VecDot(mtv, lddmm_lk, &dij2); CHKERRQ(code);

  double dij = sqrt(dij1 + dij2);

  ofstream ofile(par.dirResults() + "/D" + argv[2] + "_" + argv[3] + ".txt");
  ofile << scientific << dij;
  ofile.close();

  cout << "Distance: " << dij << endl;

  par.finalize();
  SlepcFinalize();
}
