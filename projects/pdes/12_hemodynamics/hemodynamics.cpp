/*=============================================================================
  This file is part of the code MAD
  Copyright (C) 2021-2023,
    
     Felipe Galarce at INRIA/WIAS/PUCV 

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

//#include <plot.hpp>
//#include <carotid.hpp>
//#include <phantom.hpp>
#include <toy.hpp>

int main(int argc, char *argv[]){

  Parameters par = MADinitialize(argc, argv);

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);
    
  Geometry geo;
  geo.initialize(par, io);

  Boundary bd;
  bd.initialize(par, geo);

  NavierStokes ns;
  ns.initialize(par, geo, bd);

  BoundaryCondition bc;
  bc.initialize(bd, geo, par);

//  Ploter ploter;
//  ploter.initialize(bd, geo, par);

  int nbDofs = io.nbVertices()*(par.nbDofsPerNode()[0] + par.nbDofsPerNode()[1]);

  /* initial condition */
  if (par.initial_condition_pre() != "none" && par.initial_condition_vel() != "none"){
    Vec u0 = zeros(nbDofs);
    io.loadVector(u0, par.initial_condition_vel(), par.initial_condition_pre());
    io.writeState(u0, 0.0);
    ns.setInitialCondition(u0);
  }
  Vec u = zeros(nbDofs);

  Mat A = ns.assembleLHS_static();
  ns.setSolver();

  Vec b = ns.assembleRHS();
  double TOL = 1e-1;
//  ofstream bfFile(par.results() + "/bf.txt");
//  for (int i = 0; i < par.outlets().size(); i++){
//    bfFile << par.outlets()[i] << " ";
//  }
//  bfFile << endl;
  for (double t : range(0, par.nbIterations()*par.timeStep(), par.timeStep())){

    MADprint("\n- - - - - - - - - - - - - - - - - - - - - - - - -\n");
    MADprint("Navier Stokes: Solving Navier-Stokes equation for time: ", t);

    Mat C = ns.assembleLHS(A);
    Vec b = ns.assembleRHS();
    bc.applyBC(u, C, b, t);
    ns.setLHS(C);
    ns.solve(b, u);
    bd.backflow(u);
    ns.update(u, t);
    MatDestroy(&C);
    VecDestroy(&b);

    if (t/par.timeStep() >= par.start() && ((int)floor(t/par.timeStep())) % par.jump() == 0){
      Vec velocity = ns.calculus.split(u)[0];
      io.writeState(velocity , "Velocity", t + par.timeStep());
      Vec pressure = ns.calculus.split(u)[1];
      io.writeState(pressure, "Pressure", t + par.timeStep());
//      vector<Vec> uCoord = ns.calculus.decompose_vector_field(velocity);
//      io.writeState(uCoord[0], "VelocityX", t + par.timeStep());
//      io.writeState(uCoord[1], "VelocityY", t + par.timeStep());
//      io.writeState(uCoord[2], "VelocityZ", t + par.timeStep());
//      Vec gradP = ns.calculus.gradient(pressure);
//      io.writeState(gradP, "grad(P)", t + par.timeStep());
//      Vec gradUx = ns.calculus.gradient(uCoord[0]);
//      Vec gradUy = ns.calculus.gradient(uCoord[1]);
//      Vec gradUz = ns.calculus.gradient(uCoord[2]);
//      io.writeState(gradUx, "grad(Ux)", t + par.timeStep());
//      io.writeState(gradUy, "grad(Uy)", t + par.timeStep());
//      io.writeState(gradUz, "grad(Uz)", t + par.timeStep());
    }
    ns.computeFlows(u);

//    ploter.plot(u);
  }
//bfFile.close();

  MatDestroy(&A);
  VecDestroy(&u);
  ns.finalize();
  MADfinalize(par);
}
