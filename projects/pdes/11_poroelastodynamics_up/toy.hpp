#include <mad.hpp>

vector<double> csf(vector<double> x, double t, Parameters par){
  vector<double> bc(1, par.csf_pressure());
  return bc;
}

vector<double> top(vector<double> x, double t, Parameters par){
  vector<double> bc(2, 0.0);
  bc[1] = par.amplitude();
  return bc;
}

vector<double> bottom(vector<double> x, double t, Parameters par){
  vector<double> bc(2, 0.0);
  return bc;
}

void applyBC(Boundary bd, Parameters par, Mat A){
  bd.Dirichlet(1, bottom);
  /* Uniform pressure at boundaries  */
  for (int i : par.walls()){
    bd.Dirichlet(i, csf, 1);
  }
  bd.block(A);
}

void applyBC(Boundary bd, Parameters par, Vec b){
  bd.NeumannNormalConstant(3, par.amplitude());
  bd.Dirichlet(1, bottom);
  /* Uniform pressure at boundaries  */
  for (int i : par.walls()){
    bd.Dirichlet(i, csf, 1);
  }
  bd.block(b);
}
