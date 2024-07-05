#include<ultra-4d-flow.hpp>

int main(int argc, char *argv[], char *envp[]){

  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL);

  /* Parse parameters file */
  assert(argc == 2);
  Parameters par(argv[1]);

  par.print();

  CFD cfd;
  cfd.initialize(par);

  ofstream p_drop(par.dirResults() + "/p_drop.txt");

  p_drop << scientific << 0.0 << " " << 0.0 << endl;

  for (double time : range(par.UStimeStart(), par.UStimeEnd() - par.UStimeStep(), par.UStimeStep())){

    int iteration;
    if ( abs(time - floor(time)) > 0.5){
      iteration = floor(time/0.004) + 1;
    } else {
      iteration = floor(time/0.004);
    }

    cout << "\n- - - - - - - - - - - - - -" << endl;
    cout << "Imrp iteration: " << iteration << endl;

    vector<double> pDrop = cfd.pressureDropVirtualWork(par.dirResults() + "/v_star." + wildcard(iteration+1) + ".vct", par.dirResults() + "/v_star." + wildcard(iteration) + ".vct");
  
    cout << "Pressure drop virtual work: " << endl;
    cout << "2-3: " << scientific << pDrop[0] * 760.0 / ( 10.0 * 101325.0) << endl;
    cout << "2-4: " << scientific << pDrop[1] * 760.0 / ( 10.0 * 101325.0)  << endl;

    p_drop << scientific << pDrop[0] << " " << pDrop[1] << endl;

  } 

  p_drop.close();

  cfd.finalize();
  SlepcFinalize();

}
