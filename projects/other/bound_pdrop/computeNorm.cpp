#include <ultra-4d-flow.hpp>

int main(int argc, char *argv[]){

  assert(argc == 2);

  string data_file = argv[1];
  
  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  Parameters par(data_file);
  par.print();

  /* Initialize MDMA */
  CFD cfd;
  cfd.initialize(par);
  InnerProduct ip;

  int nbVertices = cfd.io.nbVertices();
  int nbDofs = 3 * cfd.io.nbVertices(); 

  /* This are matrices for the potential problem, therefore they are talking about scalar magnitudes */
  Mat mass = mat(nbVertices, nbVertices);
  loadMat(mass, "../getMatrices/results/mass.bin");

  ofstream div_file("./norm_div_u.txt");
  for (double time : range(par.UStimeStart(), par.UStimeEnd(), par.UStimeStep())){

    int iteration;
    if ( abs(time - floor(time)) > 0.5){
      iteration = floor(time/0.004) + 1;
    } else {
      iteration = floor(time/0.004);
    }

    cout << "\n= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =" << endl;
    cout << "Iteration: " << iteration << endl;
    cout << "Computing divergence." << endl;

    Vec u = vec(nbVertices);
    loadVec(u, "./results/div_u." + wildcard(iteration) + ".scl");

    double norm_divu = sqrt(ip(u, mass, u));
    cout << "|| div u || = " << norm_divu << endl;

    div_file << norm_divu << endl;

  }
  div_file.close();
  cfd.finalize();
  SlepcFinalize();

}
