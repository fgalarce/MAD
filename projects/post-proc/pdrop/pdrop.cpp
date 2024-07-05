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

  ofstream pdrop_file("../../data/reconstruction/up_H1/sim" + wildcard(patientId) + "/pdrop.txt");
  pdrop_file << "pdrop_o1\t";
  pdrop_file << "pdrop_o2\t";
  pdrop_file << "pdrop_o1_gt\t";
  pdrop_file << "pdrop_o2_gt" << endl;

  for (int i = 0; i < nbIter; i++){

    cout << "\n- - - - - - - - - - - - - -" << endl;
    cout << "Iteration " << i << endl;

    /* Load p* */
    Vec p_star = vec(cfd.io.nbVertices());
    loadVec(p_star, "../../data/reconstruction/up_H1/sim" + wildcard(patientId) + "/p_star." + wildcard(i + itInit) + ".scl");

    /* Load p */
    Vec p = vec(cfd.io.nbVertices());
    loadVec(p, "../../data/targets/sim" + wildcard(patientId) + "/pressure." + wildcard(i + itInit) + ".scl");

    double p_inlet = cfd.integralSurface(p_star, 2) / cfd.geo.computeBoundaryArea(2);
    double p_out3 = cfd.integralSurface(p_star, 3) / cfd.geo.computeBoundaryArea(3);
    double p_out4 = cfd.integralSurface(p_star, 4) / cfd.geo.computeBoundaryArea(4);

    double p_inlet_gt= cfd.integralSurface(p, 2) / cfd.geo.computeBoundaryArea(2);
    double p_out3_gt = cfd.integralSurface(p, 3) / cfd.geo.computeBoundaryArea(3);
    double p_out4_gt = cfd.integralSurface(p, 4) / cfd.geo.computeBoundaryArea(4);

    cout << "Pressure drop: " << endl;
    cout << "2-3: " << scientific << (- p_inlet + p_out3) * 760.0 / ( 10.0 * 101325.0) << endl;
    cout << "2-4: " << scientific << (- p_inlet + p_out4) * 760.0 / ( 10.0 * 101325.0) << endl;

    pdrop_file << scientific << p_out3 - p_inlet << " ";
    pdrop_file << scientific << p_out4 - p_inlet << " ";
    pdrop_file << scientific << p_out4_gt - p_inlet_gt << " ";
    pdrop_file << scientific << p_out3_gt - p_inlet_gt << endl;

    VecDestroy(&p_star);
    VecDestroy(&p);

  } 

  pdrop_file.close();
  cfd.finalize();
  SlepcFinalize();

}
