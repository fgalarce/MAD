#include <ultra-4d-flow.hpp>

int main(int argc, char *argv[]){

  assert(argc == 2);

  string data_file = argv[1];
  
  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  Parameters par(data_file);
  par.print();

  Patient constance; constance.initialize(par);

  InnerProduct ip = constance.ip;

  double max_signal = 0;

  for (double time : range(par.UStimeStart(), par.UStimeEnd(), par.UStimeStep())){

    cout << "\n- - - - - - - - - - - - - - - " << endl;
    cout << "Time: " << time << endl;

    constance.measure(time);

    PetscReal max_signal0;
    VecMax(constance.usMeasures().measures(), NULL, &max_signal0);

    if (max_signal0 > max_signal)
      max_signal = max_signal0;

  }

  ofstream max_l(par.dirResults() + "/delta.txt");
  max_l << max_signal << endl;
  max_l.close();
  

  constance.discharge();
  SlepcFinalize();

}
