vector<double> inlet(vector<double> x, double t, Parameters par){
  vector<double> bc(2, 0.0);
//  vector<double> center(2, 0.0);
//  center[0] = -10.0;
//  center[1] = 5.0;
//  double r = 4.0;
//  bc[0] = par.amplitude()/2.0*(1.0 - dot(center - x, center - x) / (r*r) ) * (sin(2*PI * t / par.period()) + 1);
//  bc[0] = 0.001/2.0 * (cos(2*PI * t / par.period()) + 1);

  vector<double> center(2,0.0);
  double r = 3.999999999999999;
  bc[0] = -0.001 * (dot(center - x, center - x)/(r*r) - 1.0);
  return bc;
}

vector<double> inlet_vel(vector<double> x, double t, Parameters par){
  vector<double> bc(2, 0.0);
  bc[0] = 0.0;
  return bc;
}

vector<double> outlet(vector<double> x, double t, Parameters par){
  vector<double> bc(2, 0.0);
  return bc;
}

vector<double> zero(vector<double> x, double t, Parameters par){
  vector<double> bc(1, 0.0);
  return bc;
}

vector<double> wall(vector<double> x, double t, Parameters par){
  vector<double> bc(2, 0.0);
  return bc;
}
