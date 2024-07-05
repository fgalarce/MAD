#include<ultra-4d-flow.hpp>


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

  cout << "Boundary labels: " << (int)par.bdLabels(0) << endl;
  cout << "Boundary labels: " << (int)par.bdLabels(1) << endl;
  cout << "Boundary labels: " << (int)par.bdLabels(2) << endl;

  ofstream flowFile(par.dirResults() + "flow.txt");

  for (int i = 0; i < nbIter; i++){

    cout << "\n- - - - - - - - - - - - - -" << endl;
    cout << "Iteration " << i << endl;

    Vec u = vec(cfd.io.nbVertices()*3);
    loadVec(u, "../../data/targets/sim" + wildcard(patientId) + "/velocity." + wildcard(i + itInit) + ".vct");

    for (int j : par.bdLabels()){
      double flow_u = cfd.flow(u,j);
      cout << "Flow: " << flow_u << endl; 
      flowFile << scientific << flow_u << " ";
    }
    flowFile << endl;

  } 
  flowFile.close();
  cfd.finalize();
  SlepcFinalize();

}
