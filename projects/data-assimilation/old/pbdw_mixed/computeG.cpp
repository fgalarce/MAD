#include<ultra-4d-flow.hpp>

int main(int argc, char *argv[], char *envp[]){

  /* Parse parameters file */
  assert(argc == 2);
  Parameters par(argv[1]);
  par.print();

  /* Deploy petsc */
  SlepcInitialize(&argc, &argv, NULL, NULL);

  IO io;
  io.initialize(par);

  Geometry geo;
  geo.initialize(par, io);

  InnerProduct ip;
  ip.initialize(par, geo);

  int nbDofs = io.nbVertices()*4;
  int nbVertices = io.nbVertices();

  /* Import Riesz representers */
  int m = 216;
  vector<Vec> rr(m);
  for (int i = 0; i < m; i++){
    vec(rr[i], 3*nbVertices);
    loadVec(rr[i], "/local/fgalarce/4d-flow-ultrasound/data/reconstruction/rr_H1/rieszRepresenter." + wildcard(i) + ".bin");
  }

  /* Import basis */
  vector<Vec> phi(par.nbModes());
  for (int i = 0; i < par.nbModes(); i++){
    Vec vel = vec(3*nbVertices);
    Vec press = vec(nbVertices);
    loadVec(vel, par.dirResults()+ "mode_velocity." + wildcard(i) + ".vct");
    loadVec(press, par.dirResults() + "mode_pressure." + wildcard(i) + ".scl");

    vector<double> vel_arr(3*nbVertices);
    vector<double> press_arr(nbVertices);

    vec(phi[i], 4*nbVertices);
    VecGetValues(vel, 3*nbVertices, &range(3*nbVertices)[0], &vel_arr[0]);
    VecGetValues(press, nbVertices, &range(nbVertices)[0], &press_arr[0]);
    VecSetValues(phi[i], 3*nbVertices, &range(3*nbVertices)[0], &vel_arr[0], INSERT_VALUES);
    VecSetValues(phi[i], nbVertices, &range(3*nbVertices, 4*nbVertices)[0], &press_arr[0], INSERT_VALUES);
    VecAssemblyBegin(phi[i]);
    VecAssemblyEnd(phi[i]);
    
    VecDestroy(&vel);
    VecDestroy(&press);
  }

  /* Assemble G */
  Tic tic;
  for (int n : range(10,60,10)){
    cout << "n: " << n << endl;
    Mat G = mat(m, n);
    for (int i = 0; i < m; i++){
      cout << "Assemblig G: row " << i+1 << "/" << m << "\r";
      for (int j = 0; j < n; j++){
        MatSetValue(G, i, j, ip(rr[i], phi[j]), INSERT_VALUES);
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

  for (int i = 0; i < par.nbModes(); i++)
    VecDestroy(&phi[i]);  

  for (int i = 0; i < m; i++)
    VecDestroy(&rr[i]);
  SlepcFinalize();
}
