/*=============================================================================
  This file is part of the code M.D.M.A. (or MDMA)
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2021,
    
     Felipe Galarce at INRIA

  MDMA is free software; you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by the Free
  Software Foundation; either version 2.1 of the License, or (at your option)
  any later version.

  MDMA is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with MDMA. If not, see http://www.gnu.org/licenses/.
  =============================================================================*/

#include <ultra-4d-flow.hpp>

/* function pointers are passed as argument to BC class */
vector<double> inlet(vector<double> x, double t){
  vector<double> bc(3, 0.0);
//  double r = 1.0;
//  bc[1] = (1 - (x[0]*x[0] + x[1]*x[1])/(r*r));
  bc[1] = 1.0;
  return bc;
}
vector<double> noslip(vector<double> x, double t){
  vector<double> bc(3,0.0);
  return bc;
}

int main(int argc, char *argv[]){
  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  /* Parse data file */
  assert(argc == 2);
  string data_file = argv[1];
  Parameters par(data_file);
  par.print();

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);
    
  Geometry geo;
  geo.initialize(par, io);

  MasterTetrahedron tet; 
  tet.initialize(par);

  BoundaryConditions bc;
  bc.initialize(par, geo);

  int nbDofs = 4*io.nbVertices();
  int nbDofsPress = io.nbVertices();
  int nbDofsVel = 3*io.nbVertices();
  int nbVertices = io.nbVertices();

  Mat M = mat(nbDofs, nbDofs); /* Mass matrix */
  Mat B = mat(nbDofs, nbDofs); /* Static matrix */

  int m,n;
  code = MatGetOwnershipRange(M, &m, &n); CHKERRQ(code);

  Tic tic;
  for (int tetraId = 0; tetraId < geo.tetrahedron()[0].size(); tetraId++){ /* loop on tetra */

    vector<int> tetra = geo.tetrahedron()[0][tetraId];
    /* get finite element coordinates */
    vector<vector<double>> coordinates(4);
    coordinates[0] = geo.coordinates()[tetra[0]];
    coordinates[1] = geo.coordinates()[tetra[1]];
    coordinates[2] = geo.coordinates()[tetra[2]];
    coordinates[3] = geo.coordinates()[tetra[3]];

    tet.setCoordinates(coordinates);
    tet.computeSize();

    if (tetraId % (geo.tetrahedron()[0].size() / 10) == 0){
      if (world_rank == 0) cout << "  \n* * * * * * * * * * * * * * * " << endl;
      if (world_rank == 0) cout << "  Elementary matrix for tetra: " << tetraId << "/" << geo.tetrahedron()[0].size() - 1 << endl;
      if (world_rank == 0) cout << "  Coordinates: " << endl << coordinates;
      if (world_rank == 0) cout << "  Labels: " << tetra << endl; 
      if (world_rank == 0) cout << "  Size: " << tet.size() << endl;
      if (world_rank == 0) cout << "  det(J) = " << tet.detJacobian() << endl; 
      if (world_rank == 0) cout << "  volume = " << tet.volume() << " cms^3." << endl;
      if (world_rank == 0) cout << "  Mass(1,1) = " << tet.mass(1,1) << endl;
      if (world_rank == 0) cout << "  Stiff(1,2) = " << tet.stiffness(1,2) << endl;
      if (world_rank == 0) cout << "  Mixed(2,2) = " << tet.mixed(2,2) << endl;
    }

    /* Assemble elementary matrices */
    for (int j = 0; j < 4; j++){
      for (int i = 0; i < 4; i++){
        for (int comp = 0; comp < 3; comp++){
          if (3*tetra[i]+comp >= m && 3*tetra[i]+comp < n){
            /* M */
            double Mij = par.density()/par.timeStep()*tet.mass(i,j);
            /* stiffness */
            double Kij = tet.stiffness(i,j) ;
            code = MatSetValue(M, 3*tetra[i]+comp, 3*tetra[j]+comp, Mij, ADD_VALUES); CHKERRQ(code);
            code = MatSetValue(B, 3*tetra[i]+comp, 3*tetra[j]+comp, par.viscosity() * Kij, ADD_VALUES); CHKERRQ(code);
            /* div_phi_i phi_j */
            double Bij = tet.mixed(i,j,comp);
            code = MatSetValue(B, 3*tetra[i]+comp, nbDofsVel + tetra[j], -1.0*Bij, ADD_VALUES); CHKERRQ(code); /* -B*/
          } 
          if (nbDofsVel + tetra[j] >= m && nbDofsVel + tetra[j] < n){
            /* div_phi_i phi_j */
            double Bij = tet.mixed(i,j,comp);
            code = MatSetValue(B, nbDofsVel + tetra[j], 3*tetra[i]+comp, 1.0*Bij, ADD_VALUES); CHKERRQ(code); /* B^T */
          }
        }
        if (nbDofsVel + tetra[i] >= m && nbDofsVel + tetra[i] < n){
          double Kij = tet.stiffness(i,j);
          /* Brezzi-Pitkaranta stabilization for P1-P1 elements */
          code = MatSetValue(B, nbDofsVel + tetra[i], nbDofsVel + tetra[j], tet.size()*tet.size() * Kij, ADD_VALUES); CHKERRQ(code);
        }
      }
    }
  }
  if (world_rank == 0) cout << endl << "Time elapsed: " << tic.toc() << " sec. "<< endl;

  code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

  Vec phi0 = zeros(nbDofs);
  Vec phi = zeros(nbDofs);

  /* LHS static matrix */
  Mat A = mat(nbDofs, nbDofs);
  code = MatDuplicate(M, MAT_COPY_VALUES, &A); CHKERRQ(code);
  code = MatAXPY(A, 1.0, B, DIFFERENT_NONZERO_PATTERN); CHKERRQ(code);

  saveMat(B, "B.m");
  exit(1);


  /* Krylov solver */
  KSP ksp;
  code = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERR(code); 
  PetscBool reuse_prec = PETSC_TRUE;

  /* Save initial condition */
  io.writeState(phi, "phi", 0.0);

  for (double t : range(0, par.nbIterations()*par.timeStep(), par.timeStep())){

    if (world_rank == 0) cout << "\n - - - - - - - - - - - - - - - - - - - " << endl;
    if (world_rank == 0) cout << "Solving Navier-Stokes equation for time: " << t << endl;
    if (world_rank == 0) cout << "Iteration: " << t/par.timeStep() << endl;

    /* block and assembly rhs */
    Vec rhs = vec(nbDofs);
    MatMult(M, phi0, rhs);

    /* Assemble time dependent matrix */
    Mat C = mat(nbDofs, nbDofs);

    Vec phi0_seq;
    Vec phi0_copy = zeros(nbDofs);
    code = VecAXPY(phi0_copy, 1.0, phi0);
//    code = VecCopy(phi0, phi0_copy); CHKERRQ(code);

    code = VecCreateSeq(PETSC_COMM_SELF, nbDofs, &phi0_seq); CHKERRQ(code);
    code = VecSetFromOptions(phi0_seq); CHKERRQ(code); 
    VecScatter inctx;
    code = VecScatterCreateToAll(phi0_copy, &inctx, NULL); CHKERRQ(code);
    code = VecScatterBegin(inctx, phi0_copy, phi0_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(code);
    code = VecScatterEnd(inctx, phi0_copy, phi0_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(code);
    code = VecScatterDestroy(&inctx); CHKERRQ(code);
    code = VecDestroy(&phi0_copy); CHKERRQ(code);

    Tic tic;
    for (int tetraId = 0; tetraId < geo.tetrahedron()[0].size(); tetraId++){ /* loop on tetra */

      vector<int> tetra = geo.tetrahedron()[0][tetraId];
      /* get finite element coordinates */
      vector<vector<double>> coordinates(4);
      coordinates[0] = geo.coordinates()[tetra[0]];
      coordinates[1] = geo.coordinates()[tetra[1]];
      coordinates[2] = geo.coordinates()[tetra[2]];
      coordinates[3] = geo.coordinates()[tetra[3]];

      tet.setCoordinates(coordinates);
      tet.computeSize();

      if (tetraId % (geo.tetrahedron()[0].size() / 5) == 0){
        if (world_rank == 0) cout << "  \n - - - - - - - - - - - - - - - - - - - " << endl;
        if (world_rank == 0) cout << "  Elementary matrix for tetra: " << tetraId << "/" << geo.tetrahedron()[0].size() - 1 << endl;
        if (world_rank == 0) cout << "  Coordinates: " << endl << coordinates;
      }

      /* loc2glob dof mapping */
      vector<int> loc2glob(12);
      for (int i = 0; i < 4; i++){
        for (int comp = 0; comp < 3; comp++){
          loc2glob[3*i+comp] = 3*tetra[i]+comp;
        }
      }

      /* Gather advection field on element */
      vector<double> u_el(12);
      code = VecGetValues(phi0_seq, 12, &loc2glob[0], &u_el[0]); CHKERR(code);

      /* Set values from local to global */
      for (int j = 0; j < 4; j++){
        vector<double> u_node_j(3);
        for (int i = 0; i < 4; i++){
          u_node_j[0] = u_el[3*j];
          u_node_j[1] = u_el[3*j+1];
          u_node_j[2] = u_el[3*j+2];
          for (int comp = 0; comp < 3; comp++){
            if (3*tetra[i]+comp >= m && 3*tetra[i]+comp < n){
              code = MatSetValue(C, 3*tetra[i]+comp, 3*tetra[j]+comp, tet.advection(i, j, u_node_j), ADD_VALUES); CHKERRQ(code);
            }
          }
        }
      }
    }
    if (world_rank == 0) cout << endl << "Time elapsed: " << tic.toc() << " sec. "<< endl;

    code = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
    code = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

    /* Add static part */
    MatAXPY(C, 1.0, A, DIFFERENT_NONZERO_PATTERN);

    /* block and assembly C */
    bc.time(t);
//    bc.DirichletNormal(par.inlet(), par.inlet_u0());
    bc.Dirichlet(par.inlet(), inlet);
    for (int i = 0; i < par.walls().size(); i++){
      bc.Dirichlet(par.walls()[i], noslip);
    }
    bc.block(rhs);
    bc.block(C);

    code = KSPSetOperators(ksp, C, C); CHKERR(code); 
    code = KSPSetType(ksp, KSPGMRES); CHKERR(code); 
    PC pc;
    code = KSPGetPC(ksp, &pc); CHKERR(code); 
    code = PCSetType(pc, PCASM); CHKERR(code); 
    code = KSPMonitorSet(ksp, krylovMonitor, NULL, NULL); CHKERR(code);
    code = KSPSetReusePreconditioner(ksp, reuse_prec); CHKERR(code);
    code = KSPSetFromOptions(ksp); CHKERR(code); 

    if (world_rank == 0) cout << "  Solving linear system" << endl;
    code = KSPSolve(ksp, rhs, phi); CHKERRQ(code); 

    io.writeState(phi, "phi", t + par.timeStep());

    VecDestroy(&phi0_seq);

    VecZeroEntries(phi0);
    code = VecAXPY(phi0, 1.0, phi); CHKERRQ(code);
//    code = VecCopy(phi, phi0); CHKERRQ(code);

    double norm_sol = norm(phi);
    if (world_rank == 0) cout << "  Norm solution = " << norm_sol << endl;

    VecDestroy(&rhs);
    MatDestroy(&C);
  }

  KSPDestroy(&ksp);
  MatDestroy(&A);
  VecDestroy(&phi);
  VecDestroy(&phi0);
  MatDestroy(&B);
  MatDestroy(&M);
  SlepcFinalize();
}
