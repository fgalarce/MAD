#include <iostream>
#include <math.h>
#include <vector>
#include <slepc.h>
#include <STLvectorUtils.hpp>
#include <slepcWrapper.hpp>
#include <tools.hpp>
#include <parameters.hpp>
#include <cfd.hpp>
#include <assert.h>
#include <innerProduct.hpp>

using namespace std;

double normL2(Vec u){
  PetscInt size;
  VecGetSize(u, &size);
  Mat m = mat(size,size); 
  loadMat(m, "../../data/massBD.bin"); 

  double norm;
  Vec mtv = vec(size);
  MatMult(m, u, mtv);
  VecDot(mtv, u, &norm);
  
  MatDestroy(&m);
  VecDestroy(&mtv);
  return sqrt(norm);
}

int main(int argc, char *argv[], char *envp[]){

  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL);

  /* Parse parameters file */
  assert(argc == 3);
  Parameters par(argv[1]);
  int patientId = stoi(argv[2]);
  int itInit = par.UStimeStart()/par.UStimeStep();
  int nbIter = par.UStimeEnd()/par.UStimeStep() - itInit;
  par.print();

  CFD cfd;
  cfd.initialize(par);

  IO io = cfd.io();

  /* CI file */
  ofstream CI_file("./ci/wss.ci");

  ofstream error_file(par.dirResults() + "/error.txt");
  error_file << "errorL2\tnormGT\tnormL2_HE\tnormL2_GT_HE\tnormL2_average\n";

  cout << "====================================\n";
  cout << "Patient: " << patientId << endl;
  cout << "itInit: " << itInit << endl;
  cout << "nbIter: " << nbIter << endl;
  cout << "====================================\n";

  PetscInt indexes_average[io.nbVertices()*3];
  PetscReal array_average[io.nbVertices()*3];

  for (int i = 0; i < nbIter; i++){

    cout << "\n- - - - - - - - - - - - - -" << endl;
    cout << "Iteration " << i << endl;

    /* Load velocity */
    Vec u = vec(io.nbVertices()*3);
    loadVec(u, "../../data/targets/sim" + wildcard(patientId) + "/velocity." + wildcard(i + itInit) + ".vct");

    Vec v = vec(io.nbVertices()*3);
    loadVec(v, "../../data/reconstruction/u_H1/sim" + wildcard(patientId) + "/v_star." + wildcard(i) + ".vct");

    /* Compute WSS */
    Vec wss_gt = cfd.wallShearStress(u);
    Vec error = cfd.wallShearStress(v);

    /* Write to ensight */
    io.writeState(wss_gt, "wss_gt", i);
    io.writeState(error, "wss", i);

    VecAXPY(error, -1.0, wss_gt);

    Vec harmonic_error = vec(io.nbVertices()*3); 
    harmonic_error = cfd.harmonicExtensionNeumann(error);

    Vec harmonic_gt = vec(io.nbVertices()*3);
    harmonic_gt = cfd.harmonicExtensionNeumann(wss_gt);

    /* Write to ensight */
    io.writeState(harmonic_error, "HE_error", i);
    io.writeState(harmonic_gt, "HE_gt", i);

    double errorHE = normL2(cfd.trace(harmonic_error)); 
    double normGT_HT = normL2(cfd.trace(harmonic_gt));

    double errorL2 = normL2(error); 
    double normGT = normL2(wss_gt);

    /* Compute mean vector field over boundary */
    vector<double> intGonBoundary = cfd.computeVectorIntegralOverSurface(error);
    Vec average = vec(io.nbVertices()*3);
    for (int i = 0; i <  io.nbVertices(); i++){
      for (int comp = 0; comp < 3; comp++){
        indexes_average[3*i+comp] = 3*i+comp;
        array_average[3*i+comp] = intGonBoundary[comp];
      }
    }
    VecSetValues(average, io.nbVertices()*3, indexes_average, array_average, INSERT_VALUES);
    VecAssemblyBegin(average); 
    VecAssemblyEnd(average); 
    PetscReal measureBoundary = cfd.geometry().computeBoundaryArea();
    VecScale(average, 1.0/measureBoundary);
    double normL2_average = normL2(cfd.trace(average));

    error_file << scientific << setprecision(5) << errorL2 << " " << normGT << " " << errorHE << " " << normGT_HT << " " << normL2_average << endl;
    VecDestroy(&error);
    VecDestroy(&harmonic_error);
    VecDestroy(&harmonic_gt);
    VecDestroy(&average);
    VecDestroy(&wss_gt);
    VecDestroy(&u);
    VecDestroy(&v);

    CI_file << errorHE << endl;
  } 

  error_file.close();
  CI_file.close();
  cfd.finalize();
  SlepcFinalize();

}
