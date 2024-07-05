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

#include <phantom.hpp>

int main(int argc, char *argv[]){

  Parameters par = MADinitialize(argc, argv);

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);
    
  Geometry geo;
  geo.initialize(par, io);

  Boundary boundary;
  boundary.initialize(par, geo);

  Calculus calculus;
  calculus.initialize(par, geo, boundary);

  NS ns;
  ns.initialize(par, geo, boundary);

  int nbDofs = io.nbVertices()*par.nbVariables();
//  /* initial condition */
////  if (par.initial_condition_pre() != "none" && par.initial_condition_vel() != "none"){
//    Vec u0 = zeros(nbDofs);
//    vector<string> init_cond(4);
//    init_cond[0] = "./phantom2/ux.00325.scl";
//    init_cond[1] = "./phantom2/uy.00325.scl";
//    init_cond[2] = "./phantom2/uz.00325.scl";
//    init_cond[3] = "./phantom2/p.00325.scl";
////    io.loadVector(u0, par.initial_condition_vel(), par.initial_condition_pre());
//    io.loadVector(u0, init_cond);
//    io.writeState(u0, 0.0);

  Mat A = ns.assembleLHS_static();
  ns.setSolver();

//    ns.setInitialCondition(u0);
//  }
  Vec u = zeros(nbDofs);


  Aorta fluid;
  fluid.initialize(ns, par, boundary, geo, calculus, io);

  for (double t : range(par.timeStep(), par.nbIterations()*par.timeStep(), par.timeStep())){

    MADprint("\n- - - - - - - - - - - - - - - - - - - - - - - - -\n");
    MADprint("Navier Stokes: Solving Navier-Stokes equation for time: ", t);

    /* Run non-linear solver for current time step */
    fluid.solve(A, boundary, u, t, ns);
    ns.update(u, t);
    if (t/par.timeStep() >= par.start() && ((int)floor(t/par.timeStep())) % par.jump() == 0){
      io.writeState(u, t);
    }
    ns.computeFlows(u);

    if (par.power_law_n() != 1.0){
      io.writeState(ns.viscosity_p1(), "viscosity", t);
    }

    double normSOL = norm(u);
    MADprint("Navier Stokes: norm solution = ", normSOL);

  }

  MatDestroy(&A);
  VecDestroy(&u);
  ns.finalize();
  MADfinalize(par);
}
