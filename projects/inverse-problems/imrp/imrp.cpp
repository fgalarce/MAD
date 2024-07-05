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
#include <fstream>

using namespace std;

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

  ofstream p_drop(par.dirResults() + "/p_drop.txt");
  ofstream p_dropGT(par.dirResults() + "/p_drop_GT.txt");


  p_dropGT << scientific << "0 0" << endl;
  p_drop << scientific << "0 0" << endl;

  for (int i = 1; i < nbIter; i++){

    cout << "\n- - - - - - - - - - - - - -" << endl;
    cout << "Iteration " << i << endl;

    vector<double> pDrop = cfd.pressureDropVirtualWork( "../../data/targets/sim" + wildcard(patientId) + "/velocity." + wildcard(i + itInit) + ".vct",
                                                        "../../data/targets/sim" + wildcard(patientId) + "/velocity." + wildcard(i + itInit-1) + ".vct");
  
    /* Load pressure */
    Vec p = vec(cfd.io.nbVertices());
    loadVec(p, "../../data/targets/sim" + wildcard(patientId) + "/pressure." + wildcard(i + itInit) + ".scl");

    double p_inlet = cfd.integralSurface(p, 2) / cfd.geo.computeBoundaryArea(2);
    double p_out3 = cfd.integralSurface(p, 3) / cfd.geo.computeBoundaryArea(3);
    double p_out4 = cfd.integralSurface(p, 4) / cfd.geo.computeBoundaryArea(4);

    cout << "Pressure drop: " << endl;
    cout << "2-3: " << scientific << (p_inlet - p_out3) * 760.0 / ( 10.0 * 101325.0) << endl;
    cout << "2-4: " << scientific << (p_inlet - p_out4) * 760.0 / ( 10.0 * 101325.0) << endl;

    cout << "Pressure drop virtual work: " << endl;
    cout << "2-3: " << scientific << pDrop[0] * 760.0 / ( 10.0 * 101325.0) << endl;
    cout << "2-4: " << scientific << pDrop[1] * 760.0 / ( 10.0 * 101325.0)  << endl;

    p_dropGT << scientific << p_out3 - p_inlet << " " << p_out4 - p_inlet << endl;
    p_drop << scientific << pDrop[0] << " " << pDrop[1] << endl;

  } 

  p_dropGT.close();
  p_drop.close();

  cfd.finalize();
  SlepcFinalize();

}
