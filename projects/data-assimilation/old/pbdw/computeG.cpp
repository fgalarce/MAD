#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <slepc.h>
#include <STLvectorUtils.hpp>
#include <slepcWrapper.hpp>
#include <tools.hpp>
#include <assert.h>
#include <parameters.hpp>
#include <io.hpp>

using namespace std;

int main(int argc, char *argv[], char *envp[]){

  /* Parse parameters file */
  assert(argc == 2);
  Parameters par(argv[1]);
  par.print();

  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL);

  IO io;
  io.initialize(par);
  size_t nbDofs = io.nbVertices()*3;

  /* Import Riesz representers */
  int m = 216;
  vector<Vec> rr(m);
  for (int i = 0; i < m; i++){
    vec(rr[i], nbDofs);
    loadVec(rr[i], "/local/fgalarce/4d-flow-ultrasound/data/reconstruction/rr_H1/rieszRepresenter." + wildcard(i) + ".bin");
  }

  /* Import mass and stiffnes matrix */
  Mat M = mat(nbDofs, nbDofs);
  Mat K = mat(nbDofs, nbDofs);
  loadMat(M, par.mass_matrix_path());
  loadMat(K, par.stiffness_matrix_path());
  MatAXPY(M, 1.0, K, SAME_NONZERO_PATTERN);

  /* Import basis */
  vector<Vec> phi(par.nbModes());
  for (int i = 0; i < par.nbModes(); i++){
    vec(phi[i], nbDofs);
    loadVec(phi[i], par.dirResults() + "mode." + wildcard(i) + ".vct");
  }

  /* Assemble G */
  Tic tic;
  for (int n : range(10,120,10)){
    cout << "n: " << n << endl;
    Mat G = mat(m, n);
    for (int i = 0; i < m; i++){
      cout << "Assemblig G: row " << i+1 << "/" << m << "\r";
      for (int j = 0; j < n; j++){
        double Gij;
        Vec mtv = vec(nbDofs);
        MatMult(M, phi[j], mtv);
        VecDot(rr[i], mtv, &Gij);
        VecDestroy(&mtv);
        MatSetValue(G, i, j, Gij, INSERT_VALUES);
      }
    }
    cout << endl;
    cout << "Time elapsed: " << tic.toc() << endl;

    MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY);
    saveMat(G, par.dirResults() + "/G_" + to_string(n) + "modes.bin");  

    cout << "Norm G " << norm(G) << endl;
    MatDestroy(&G);
  }


  for (int i = 0; i < m; i++)
    VecDestroy(&rr[i]);
  SlepcFinalize();
}
