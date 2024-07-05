#include<ultra-4d-flow.hpp>

int main(int argc, char *argv[]){

  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  /* Parse data file */
  assert(argc == 2);
  Parameters par(argv[1]);

  /* Initialize MAD objects */
  par.print();

  IO io;
  io.initialize(par);

  Geometry geo;
  geo.initialize(par, io);

  Learning learning;
  learning.initialize(par, geo);

  /* Dissimilarity matrix */
  Mat D = mat(par.nbMeshes(), par.nbMeshes(), "dense");
  loadMat(D, par.dissimilarity_matrix());

  /* Computing low dimensional representation of data */
  Mat mds = learning.MDS(D);

  par.finalize();
  SlepcFinalize();
}
