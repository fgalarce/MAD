#include<ultra-4d-flow.hpp>

int main(int argc, char *argv[]){

  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  /* Parse data file */
  assert(argc == 4);
//  assert(argc == 4);
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
    string field_name;
    if (stoi(argv[2]) == stoi(argv[3])){
      field_name = "/mode.";
    } else {
      field_name = "/I_PTpod.";
    }
    loadVec(phi1[i], par.templateModel() + field_name + wildcard(i) + ".vct");
  }

  if (!la.checkOrthonormality(phi1)){
    cout << "MAIN: basis from " << par.templateModel() << " are not orthonormal." << endl;
    exit(1);
  } 
  if (!la.checkOrthonormality(phi0)){
    cout << "MAIN: basis from " << par.model() << " are not orthonormal." << endl;
    exit(1);
  }

  Mat G = mat(n,n);
  for (int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      code = MatSetValue(G,i,j,ip(phi0[i], phi1[j]),INSERT_VALUES); CHKERRQ(code);
    }
  }
  code = MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  
  vector<double> sv = getSingularValues(G, n);

  cout << "MAIN: Singular values of CG: " << endl << sv << endl;
//  sv.erase(min_element(sv.begin(), sv.end()));
  double dH = sqrt(1.0 - min(sv)*min(sv));
  cout << "MAIN: Hausdorff distance = " << dH << endl;
//  cout << "MAIN: Hausdorff distance = " << 1.0 - min(sv)*min(sv)<< endl;
//  cout << "asdasd: " << sqrt(1.0 - min(sv)) << endl;

  ofstream haus(par.dirResults() + "/haus" + argv[2] + "_"+ argv[3] + ".txt");
//  ofstream haus(par.dirResults() + "/haus" + argv[2] + ".txt");
  haus << dH << endl;
  haus.close();
  
  MatDestroy(&G);
  SlepcFinalize();
}
