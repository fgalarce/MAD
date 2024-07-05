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

#include <toy.hpp>

int main(int argc, char *argv[]){

  Parameters par = MADinitialize(argc, argv);

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);
    
  Geometry geo;
  geo.initialize(par, io);

  Boundary boundary;
  boundary.initialize(par, geo);

  NavierStokes ns;
  ns.initialize(par, geo, boundary);

  BoundaryCondition bc;
  bc.initialize(boundary, geo, par);

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

  for (double t : range(0, par.nbIterations()*par.timeStep(), par.timeStep())){

    MADprint("\n- - - - - - - - - - - - - - - - - - - - - - - - -\n");
    MADprint("Navier Stokes: Solving Navier-Stokes equation for time: ", t);

    /* assemble time dependent LHS and block full matrix and RHS */
    Mat C = ns.assembleLHS(A);
    Vec b = ns.assembleRHS();
    if (par.backflowStab()){
      boundary.backflow(u);
    }
    bc.applyBC(u, C, b, t);

    /* solve linear system */
    ns.setLHS(C);
    ns.solve(b, u);
    if (par.backflowStab()){
      boundary.backflow(u);
    }
    ns.update(u, t);

    if (t/par.timeStep() >= par.start() && ((int)floor(t/par.timeStep())) % par.jump() == 0){
      io.writeState(u , t + par.timeStep());
    }
    ns.computeFlows(u);

    double normRHS = norm(b);
    double normLHS = norm(A);
    double normSOL = norm(u);
    MADprint("Navier Stokes: norm RHS blocked = ", normRHS);
    MADprint("Navier Stokes: norm LHS blocked = ", normLHS);
    MADprint("Navier Stokes: norm solution = ", normSOL);

    MatDestroy(&C);
    VecDestroy(&b);
  }

  MatDestroy(&A);
  VecDestroy(&u);
  ns.finalize();
  MADfinalize(par);
}
