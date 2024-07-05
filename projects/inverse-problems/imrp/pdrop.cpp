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

  Ensight io = cfd.io();

  vector<double> pout3;
  vector<double> pout4;

  for (int i = 0; i < nbIter; i++){

    cout << "\n- - - - - - - - - - - - - -" << endl;
    cout << "Iteration " << i << endl;

    /* Load velocity */
    Vec p = vec(io.nbVertices());
    loadVec(p, "../../data/targets/sim" + wildcard(patientId) + "/pressure." + wildcard(i + itInit) + ".scl");

    double p_inlet = cfd.integralSurface(p, 2) / cfd.geometry().computeBoundaryArea(2);
    double p_out3 = cfd.integralSurface(p, 3) / cfd.geometry().computeBoundaryArea(3);
    double p_out4 = cfd.integralSurface(p, 4) / cfd.geometry().computeBoundaryArea(4);

    cout << "Pressure drop: " << endl;
    cout << "2-3: " << scientific << abs(abs(p_inlet) - abs(p_out3)) * 760.0 / ( 10.0 * 101325.0) << endl;
    cout << "2-4: " << scientific << abs(abs(p_inlet) - abs(p_out4)) * 760.0 / ( 10.0 * 101325.0) << endl;

    pout3.push_back(p_out3 - p_inlet);
    pout4.push_back(p_out4 - p_inlet);

  } 

  exportData(pout3, par.dirResults() + "./pout3.txt");
  exportData(pout4, par.dirResults() + "./pout4.txt");

  cfd.finalize();
  SlepcFinalize();

}
