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

  ofstream kinetic(par.dirResults() + "/kinetic.txt");
  ofstream viscous(par.dirResults() + "/viscous.txt");
  ofstream viscousBD(par.dirResults() + "/viscousBD.txt");
  ofstream convective(par.dirResults() + "/convective.txt");

  for (int i = 0; i < nbIter; i++){

    cout << "\n- - - - - - - - - - - - - -" << endl;
    cout << "Iteration " << i << endl;

    double kinetic_energy = cfd.kinetic_energy("../../data/targets/sim" + wildcard(patientId) + "/velocity." + wildcard(i + itInit) + ".vct");
    cout << "Kinetic energy: " << kinetic_energy << endl; 
    kinetic << kinetic_energy << endl;

    double viscous_energy = cfd.viscous_energy("../../data/targets/sim" + wildcard(patientId) + "/velocity." + wildcard(i + itInit) + ".vct");
    cout << "Viscous energy: " << viscous_energy << endl; 
    viscous << viscous_energy << endl;  

    double convective_energy = cfd.convective_energy("../../data/targets/sim" + wildcard(patientId) + "/velocity." + wildcard(i + itInit) + ".vct");
    cout << "Convective energy: " << convective_energy << endl; 
    convective << convective_energy << endl;  

    double viscous_energyBD = cfd.viscous_energyBD("../../data/targets/sim" + wildcard(patientId) + "/velocity." + wildcard(i + itInit) + ".vct");
    cout << "Viscous energy BD: " << viscous_energyBD << endl;
    viscousBD << viscous_energyBD << endl;

  } 

  convective.close();
  kinetic.close();
  viscous.close();
  viscousBD.close();
  cfd.finalize();
  SlepcFinalize();

}
