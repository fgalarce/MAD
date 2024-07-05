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

  ofstream pdrop_file(par.dirResults() + "/pdrop.txt");
  pdrop_file << "pdrop_o1\t";
  pdrop_file << "pdrop_o2\n";

  for (double time : range(par.UStimeStart(), par.UStimeEnd(), par.UStimeStep())){

    int iteration;
    if ( abs(time - floor(time)) > 0.5){
      iteration = floor(time/0.004) + 1;
    } else {
      iteration = floor(time/0.004);
    }

    cout << "\n- - - - - - - - - - - - - -" << endl;

    /* Load p* */
    Vec p_star = vec(cfd.io.nbVertices());
    loadVec(p_star, par.dirResults() + "/p_star." + wildcard(iteration) + ".scl");

    double p_inlet = cfd.integralSurface(p_star, 2) / cfd.geo.computeBoundaryArea(2);
    double p_out3 = cfd.integralSurface(p_star, 3) / cfd.geo.computeBoundaryArea(3);
    double p_out4 = cfd.integralSurface(p_star, 4) / cfd.geo.computeBoundaryArea(4);

    cout << "Pressure drop: " << endl;
    cout << "2-3: " << scientific << (- p_inlet + p_out3) * 760.0 / ( 10.0 * 101325.0) << endl;
    cout << "2-4: " << scientific << (- p_inlet + p_out4) * 760.0 / ( 10.0 * 101325.0) << endl;

    pdrop_file << scientific << p_out3 - p_inlet << " ";
    pdrop_file << scientific << p_out4 - p_inlet << "\n";

    VecDestroy(&p_star);

  } 

  pdrop_file.close();
  cfd.finalize();
  SlepcFinalize();

}
