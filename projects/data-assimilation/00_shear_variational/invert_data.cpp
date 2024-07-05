/*=============================================================================
  This file is part of the code MAD
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2017/2021,
    
     Felipe Galarce

  MAD is free software; you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by the Free
  Software Foundation; either version 2.1 of the License, or (at your option)
  any later version.

  MAD is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with MAD. If not, see http://www.gnu.org/licenses/.
  =============================================================================*/

#include <cube.hpp>
#include <plot.hpp>

vector<double> empotrado(vector<double> x, double t, Parameters par){
  vector<double> bc(par.nbDofsPerNode()[0], 0.0);
  return bc;
}

int main(int argc, char *argv[]){

  Parameters par = MADinitialize(argc, argv);
  par.print();

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);
    
  Geometry geo;
  geo.initialize(par, io);

  InverseProblemsTools ip_tools;
  ip_tools.initialize(par, geo);

  ip_tools.assembleMassAndStiffness();

  Boundary boundary;
  boundary.initialize(par, geo);

  for (int i=0; i<par.walls().size(); i++){
    boundary.Dirichlet(par.walls()[i], empotrado);
  }
  boundary.block(ip_tools.M);
  boundary.block(ip_tools.K);

  int nbDofs = geo.nbVertices * par.nbDofsPerNode()[0];
  Vec u_test = zeros(nbDofs);

  /* Transducer objects */
  InnerProduct ip;
  ip.initialize(par, geo);

  LinearAlgebra la;
  la.initialize(par, ip);

  Measures measures;
  measures.initialize(par, geo, la, io);

//  io.loadVector(u_test, par.dirModel());
  ofstream fileShear;
  fileShear.open("./G.txt");

  /* compute quadratic forms and G */
  for (double t : range(2*par.timeStep(), par.nbIterations()*par.timeStep(), par.timeStep())){

    int iteration = t/par.timeStep();
    MADprint("\n - - - - - - - - - - - - - - - - - - - \n");
    MADprint("Solving inverse problem for time: ", t);
    MADprint("Iteration: ", iteration);

    io.loadVector(u_test, par.dirSyntheticField() + "/" + par.variableName()[0] + "." + wildcard(iteration) + ".vct");

    /* Load "measures" */
    Vec u0 = zeros(nbDofs);
    Vec u1 = zeros(nbDofs);
    Vec u2 = zeros(nbDofs);
    io.loadVector(u0, par.dirSyntheticField() + "/" + par.variableName()[0] + "." + wildcard(iteration-2) + ".vct");
    io.loadVector(u1, par.dirSyntheticField() + "/" + par.variableName()[0] + "." + wildcard(iteration-1) + ".vct");
    io.loadVector(u2, par.dirSyntheticField() + "/" + par.variableName()[0] + "." + wildcard(iteration) + ".vct");

    double stiffTerm = quadraticForm(u_test, ip_tools.stiffness(), u2);

    Vec source = computeSource(par, geo, t);
    boundary.block(source);
    double sourceTerm = quadraticForm(u_test, ip_tools.mass(), source);
    sourceTerm = sourceTerm / stiffTerm;

    /* Compute inertia term */
    Vec ddu = zeros(nbDofs);
    VecAXPY(ddu,  1.0, u2);
    VecAXPY(ddu, -2.0, u1);
    VecAXPY(ddu,  1.0, u0);

    double inertiaTerm = quadraticForm(u_test, ip_tools.mass(), ddu);
    inertiaTerm = inertiaTerm*par.density()/(par.timeStep()*par.timeStep()*stiffTerm);

    double shearModulus = sourceTerm - inertiaTerm;
    MADprint("Inertia term: ", inertiaTerm);
    MADprint("sourceTerm: ", sourceTerm);
    MADprint("stiffTerm: ", stiffTerm);
    MADprint("Computed shear modulus: ", shearModulus);
    fileShear << shearModulus << endl;
  }
  fileShear.close();
  SlepcFinalize();
}
