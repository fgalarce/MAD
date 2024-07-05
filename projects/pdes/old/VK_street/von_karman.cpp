/*=============================================================================
  This file is part of the code M.D.M.A. (or MDMA)
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2017 - 2021,
    
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

#include <mad.hpp>

/* function pointers are passed as argument to5BC class */
vector<double> inlet(vector<double> x, double t){
  vector<double> bc(2, 0.0);
  bc[0] = 5*(1 - (x[1]*x[1])/(8*8));
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

  /* local values */
  int nbDofsPerNode = par.nbDofsPerNode()[0]; 
  int nbNodesPerElement = geo.dimension()+1; 

  /* global values */
  int nbDofs = io.nbVertices() * (nbDofsPerNode + 1);
  int nbDofsVel = io.nbVertices() * nbDofsPerNode;
  int nbDofsPress = io.nbVertices();

  MasterTriangle2D tria; 
  tria.initialize(par, nbDofsPerNode);

  BoundaryConditions bc;
  bc.initialize(par, geo);

  Vec u0 = zeros(nbDofs);
  Vec u = zeros(nbDofs);

  // WARNING This is a trick to avoid problem with boundary wrong labeling, to be corrected in cleaner way
  bc.Dirichlet(0, noslip); /* condition over obstacle */
  bc.Dirichlet(3, noslip); /* condition over obstacle */
  bc.Dirichlet(par.inlet(), inlet);
  vector<double> noslip(nbDofsPerNode, 0.0);
  for (int i = 0; i < par.walls().size(); i++){
    bc.Dirichlet(par.walls()[i], noslip);
  }
  bc.block(u0);
  io.writeState(u0, 0.0);
  Mat M = mat(nbDofs, nbDofs); /* Mass matrix */
  Mat B = mat(nbDofs, nbDofs); /* Static matrix */

  int m,n;
  code = MatGetOwnershipRange(M, &m, &n); CHKERRQ(code);

  if (world_rank == 0) cout << "NS: Assembling discretization matrix" << endl;

  for (int triaId = 0; triaId < geo.triangles()[0].size(); triaId++){ /* loop on tria */
    vector<int> triangle = geo.triangles()[0][triaId];

    /* get finite element coordinates */
    vector<vector<double>> coordinates(nbNodesPerElement);
    for (int i=0; i < nbNodesPerElement; i++){
      coordinates[i] = geo.coordinates()[triangle[i]];
    }
    tria.setCoordinates(coordinates);
    tria.computeSize();

    /* Assemble elementary matrices */
    for (int i = 0; i < nbNodesPerElement; i++){
      for (int j = 0; j < nbNodesPerElement; j++){
        for (int comp = 0; comp < nbDofsPerNode; comp++){
          if (nbDofsPerNode*triangle[i]+comp >= m && nbDofsPerNode*triangle[i]+comp < n){
            double Mij = tria.mass(i,j);
            double Kij = tria.stiffness_symmetric(i,j);
            double Bij = tria.mixed(i,j,comp);
            code = MatSetValue(M, nbDofsPerNode*triangle[i]+comp, nbDofsPerNode*triangle[j]+comp, par.density()/par.timeStep() * Mij, ADD_VALUES); CHKERRQ(code);
            code = MatSetValue(B, nbDofsPerNode*triangle[i]+comp, nbDofsPerNode*triangle[j]+comp, par.viscosity() * Kij, ADD_VALUES); CHKERRQ(code);
            code = MatSetValue(B, nbDofsPerNode*triangle[i]+comp, nbDofsVel + triangle[j], -1.0*Bij, ADD_VALUES); CHKERRQ(code); /* -B*/
          }
          if (nbDofsVel + triangle[j] >= m && nbDofsVel + triangle[j] < n){
            /* div_phi_i phi_j */
            double Bij = tria.mixed(i,j,comp);
            code = MatSetValue(B, nbDofsVel + triangle[j], nbDofsPerNode * triangle[i] + comp, 1.0*Bij, ADD_VALUES); CHKERRQ(code); /* B^T */
          }
        }
        if (nbDofsVel + triangle[i] >= m && nbDofsVel + triangle[i] < n){
          double Kij = tria.stiffness(i,j);
          /* Brezzi-Pitkaranta stabilization for P1-P1 elements */
          code = MatSetValue(B, nbDofsVel + triangle[i], nbDofsVel + triangle[j], tria.size()*tria.size() * Kij, ADD_VALUES); CHKERRQ(code);
        }
      }
    }

    if (triaId % (geo.triangles()[0].size() / 10) == 0){
      if (world_rank == 0) cout << "  \n* * * * * * * * * * * * * * * " << endl;
      if (world_rank == 0) cout << "  Elementary matrix for tria: " << triaId << "/" << geo.triangles()[0].size() - 1 << endl;
      if (world_rank == 0) cout << "  Coordinates: " << endl << coordinates;
      if (world_rank == 0) cout << "  Labels: " << triangle << endl; 
      if (world_rank == 0) cout << "  Size: " << tria.size() << endl;
      if (world_rank == 0) cout << "  Mass(1,1) = " << tria.mass(1,1) << endl;
      if (world_rank == 0) cout << "  Stiff(1,2) = " << tria.stiffness(1,2) << endl;
      if (world_rank == 0) cout << "  Mixed(2,2,1) = " << tria.mixed(2,2,1) << endl;
      if (world_rank == 0) cout << "  |Jacobian| = " << tria.detJacobian() << endl;
    }
  }

  code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

  /* LHS static matrix */
  Mat A = mat(nbDofs, nbDofs);
  code = MatDuplicate(M, MAT_COPY_VALUES, &A); CHKERRQ(code);
  code = MatAXPY(A, 1.0, B, DIFFERENT_NONZERO_PATTERN); CHKERRQ(code);

  /* Krylov solver */
  KSP ksp;
  code = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERR(code); 
  PetscBool reuse_prec = PETSC_TRUE;

  ofstream ksp_output(par.dirResults() + "/ksp.log");

  for (double t : range(0, par.nbIterations()*par.timeStep(), par.timeStep())){

    if (world_rank == 0) cout << "\n - - - - - - - - - - - - - - - - - - - " << endl;
    if (world_rank == 0) cout << "Solving NS equations for time: " << t << endl;
    if (world_rank == 0) cout << "Iteration: " << t/par.timeStep() << endl;

    /* block and assembly rhs */
    Vec rhs = vec(nbDofs);
    MatMult(M, u0, rhs);

    /* Assembling convection */
    Mat C = mat(nbDofs, nbDofs);

    Vec u0_seq;
    Vec u0_copy = zeros(nbDofs);
    code = VecAXPY(u0_copy, 1.0, u0);

    code = VecCreateSeq(PETSC_COMM_SELF, nbDofs, &u0_seq); CHKERRQ(code);
    code = VecSetFromOptions(u0_seq); CHKERRQ(code); 
    VecScatter inctx;
    code = VecScatterCreateToAll(u0_copy, &inctx, NULL); CHKERRQ(code);
    code = VecScatterBegin(inctx, u0_copy, u0_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(code);
    code = VecScatterEnd(inctx, u0_copy, u0_seq, INSERT_VALUES, SCATTER_FORWARD); CHKERRQ(code);
    code = VecScatterDestroy(&inctx); CHKERRQ(code);
    code = VecDestroy(&u0_copy); CHKERRQ(code);

    for (int triaId = 0; triaId < geo.triangles()[0].size(); triaId++){ 

      vector<int> triangle = geo.triangles()[0][triaId];
      /* get finite element coordinates */
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int i=0; i < nbNodesPerElement; i++){
        coordinates[i] = geo.coordinates()[triangle[i]];
      }
      tria.setCoordinates(coordinates);

      /* loc2glob dof mapping */
      vector<int> loc2glob(nbNodesPerElement*nbDofsPerNode);
      for (int i = 0; i < nbNodesPerElement; i++){
        for (int comp = 0; comp < nbDofsPerNode; comp++){
          loc2glob[nbDofsPerNode*i+comp] = nbDofsPerNode*triangle[i]+comp;
        }
      }

      /* Gather advection field on element */
      vector<double> u_el(nbNodesPerElement*nbDofsPerNode);
      code = VecGetValues(u0_seq, nbNodesPerElement*nbDofsPerNode, &loc2glob[0], &u_el[0]); CHKERR(code);

      if (triaId % (geo.triangles()[0].size() / 2) == 0){
        if (world_rank == 0) cout << "  \n * * * * * * * * * * * * * * * * * * * " << endl;
        if (world_rank == 0) cout << "  Elementary matrix for triangle: " << triaId << "/" << geo.triangles()[0].size() - 1 << endl;
        if (world_rank == 0) cout << "  Coordinates: " << endl << coordinates;
        if (world_rank == 0) cout << "  |Jacobian|: " << tria.detJacobian() << endl;
      }

      /* Set values from local to global */
      for (int i = 0; i < nbNodesPerElement; i++){
        for (int j = 0; j < nbNodesPerElement; j++){
          vector<double> u_node_j(nbDofsPerNode);
          for (int comp = 0; comp < nbDofsPerNode; comp++){
            u_node_j[comp] = u_el[nbDofsPerNode*j+comp];
          }
          for (int comp = 0; comp < nbDofsPerNode; comp++){
            if (nbDofsPerNode*triangle[i]+comp >= m && nbDofsPerNode*triangle[i]+comp < n){
              code = MatSetValue(C, nbDofsPerNode*triangle[i]+comp, nbDofsPerNode*triangle[j]+comp, tria.advection(i, j, u_node_j), ADD_VALUES); CHKERRQ(code);
            }
          }
        }
      }
    }
    code = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
    code = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

    /* Add static part */
    MatAXPY(C, 1.0, A, DIFFERENT_NONZERO_PATTERN);

    bc.block(rhs);
    bc.block(C);

    code = KSPSetOperators(ksp, C, C); CHKERR(code); 
    code = KSPSetType(ksp, KSPPREONLY); CHKERR(code); 
    PC pc;
    code = KSPGetPC(ksp, &pc); CHKERR(code); 
    code = PCSetType(pc, PCLU); CHKERR(code); 
    code = KSPSetReusePreconditioner(ksp, reuse_prec); CHKERR(code);
    code = KSPSetFromOptions(ksp); CHKERR(code); 

    if (world_rank == 0) cout << "NS: Solving linear system" << endl;
    code = KSPSolve(ksp, rhs, u); CHKERRQ(code); 
    int its;
    double rnorm;
    code = KSPGetResidualNorm(ksp,&rnorm); CHKERRQ(code);
    code = KSPGetIterationNumber(ksp,&its); CHKERRQ(code);
    if (world_rank == 0) cout << "KSP: " << its << " iterations. Residual norm = " << rnorm << endl;
    ksp_output << its << " " << rnorm << endl;

    if (world_rank == 0) cout << "NS: Back-flow stabilization." << endl;
    for (int i : geo.boundaryNodes(par.outlet())){
      for (int comp = 0; comp < nbDofsPerNode; comp++){
        if (nbDofsPerNode*i+comp >= m && nbDofsPerNode*i+comp < n){
          double uii;
          vector<int> index(1, nbDofsPerNode*i);
          code = VecGetValues(u, 1, &index[0], &uii);
          if (uii < 0){
            code = VecSetValue(u, nbDofsPerNode*i+comp, 0.0, INSERT_VALUES); CHKERRQ(code);
          }
        }
      }
    } 
    code = VecAssemblyBegin(u); CHKERR(code);
    code = VecAssemblyEnd(u); CHKERR(code);

    io.writeState(u, t + par.timeStep()); 

    VecZeroEntries(u0);
    code = VecAXPY(u0, 1.0, u); CHKERRQ(code);

    double norm_sol = norm(u);
    if (world_rank == 0) cout << "  Norm solution = " << norm_sol << endl;

    VecDestroy(&rhs);
    MatDestroy(&C);
  }

  ksp_output.close();

  KSPDestroy(&ksp);
  MatDestroy(&A);
  VecDestroy(&u);
  VecDestroy(&u0);
  MatDestroy(&B);
  MatDestroy(&M);
  SlepcFinalize();
  par.finalize();
}
