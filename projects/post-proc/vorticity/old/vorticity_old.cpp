/*=============================================================================
  This file is part of the code M.D.M.A. (or MDMA)
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2017/2021,
    
     Felipe Galarce at INRIA / WIAS

  MDMA is free software; you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by the Free
  Software Foundation; either version 2.1 of the License, or (at your option)
  any later version.

  MDMA is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with MDMA. If not, see http://www.gnu.org/licenses/.
  =============================================================================*/

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

double ip(Vec u, Mat M, Vec v){
  PetscInt size;
  VecGetSize(u, &size);
  PetscReal ip_val;
  Vec mtv;
  vec(mtv, size);
  MatMult(M, v, mtv);
  VecDot(u, mtv, &ip_val);
  VecDestroy(&mtv);
  return ip_val;
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

  Ensight io = cfd.io();

  /* CI file */
  ofstream CI_file("./ci/vorticity.ci");

  ofstream error_file(par.dirResults() + "/error.txt");
  error_file << "errorL2\tnormGT\n";

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

    cout << "Computing vorticity for u^*" << endl;
    Vec vorticity = cfd.vorticity(v);
    cout << "Computing vorticity for u" << endl;
    Vec vorticity_gt = cfd.vorticity(u);
    io.writeState(vorticity_gt, "vorticity_gt", i);
    io.writeState(vorticity, "vorticity", i);

    /* Compute error */
    double norm_gt = sqrt(ip(vorticity_gt, cfd.mass(), vorticity_gt));
    Vec error = vec(io.nbVertices()*3);
    VecCopy(vorticity_gt, error);
    VecAXPY(error, -1.0, vorticity);
    double norm_error = sqrt(ip(error, cfd.mass(), error));

    cout << "Relative error: " << scientific << setprecision(5) << norm_error / norm_gt << endl;
    error_file << scientific << setprecision(5) << norm_error << " " << norm_gt << endl;

    VecDestroy(&vorticity_gt);
    VecDestroy(&vorticity);
    VecDestroy(&u);
    VecDestroy(&v);
  } 

  error_file.close();

  cfd.finalize();
  SlepcFinalize();

}
