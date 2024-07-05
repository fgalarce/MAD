vector<double> inlet(vector<double> x, double t, Parameters par){
//  bc[0] = par.amplitude();

  vector<double> bc(3, 0.0);
  vector<double> center(3,0.0);
  center[1] = 3.0;
  center[2] = 3.0;
  double r = 3.0;
  bc[0] = -par.amplitude() * (dot(center - x, center - x)/(r*r) - 1.0);
  if (bc[0] < 0){
    bc[0] = 0.0;
  }
  return bc;
}

vector<double> outlet(vector<double> x, double t, Parameters par){
  vector<double> bc(3, 0.0);
  return bc;
}
