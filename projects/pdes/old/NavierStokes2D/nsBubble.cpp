/*=============================================================================
  This file is part of the code MAD
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2017 - 2021,
    
     Felipe Galarce at INRIA

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

#include <mad.hpp>

/* function pointers are passed as argument to BC class */
vector<double> inlet(vector<double> x, double t, Parameters par){
  vector<double> bc(2, 0.0);
  vector<double> center(2, 0.0);
  center[1] = 20.0;
  double r = 20;
  bc[0] = (1.0 - dot(center - x, center - x) / (r*r));
  return bc;
}

vector<double> noslip(vector<double> x, double t, Parameters par){
  vector<double> bc(2, 0.0);
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
  int nbVertexPerElement = geo.dimension()+1; 
  int nbNodesPerElementVelocity = geo.dimension() + 2;
  int nbNodesPerElementPressure = geo.dimension() + 1;

  /* global values */
  int nbBubbles = io.nbVertices();
  int nbDofsVelP1 = io.nbVertices() * nbDofsPerNode;
  int nbDofsVel = nbDofsVelP1 + nbBubbles;
  int nbDofsPress = io.nbVertices();
  int nbDofs = nbDofsVel + nbDofsPress;

  MasterTriangle2D_mixed tria; 
  tria.initialize(par);

  Boundary boundary;
  boundary.initialize(par, geo);

  Vec u0 = zeros(nbDofs);
  Vec u = zeros(nbDofs);
  boundary.Dirichlet(par.inlet(), inlet, 0);
  for (int i = 0; i < par.walls().size(); i++){
    boundary.Dirichlet(par.walls()[i], noslip, 0);
  }
  boundary.block(u0);
//  io.writeState(u0, 0.0);
  Mat M = mat(nbDofs, nbDofs); /* Mass matrix */
  Mat B = mat(nbDofs, nbDofs); /* Static matrix */

  int m,n;
  code = MatGetOwnershipRange(M, &m, &n); CHKERRQ(code);

  if (world_rank == 0) cout << "NS: Assembling discretization matrix" << endl;
  /*
  | A    -B |
  | B^t   0 | 
  */

  /* Some petsc solver requires explicit entries in the diagonal, so this is necessary for saddle point problems */
  MatDiagonalSet(M, zeros(nbDofs), INSERT_VALUES);
  MatDiagonalSet(B, zeros(nbDofs), INSERT_VALUES);

  for (int triaId = 0; triaId < geo.triangles()[0].size(); triaId++){ /* loop on tria */
    vector<int> triangle = geo.triangles()[0][triaId];

    /* gather coordinates and compute Jacobian */
    vector<vector<double>> coordinates(nbVertexPerElement);
    for (int i=0; i < nbVertexPerElement; i++){
      coordinates[i] = geo.coordinates()[triangle[i]];
    }
    tria.setCoordinates(coordinates);
    tria.computeSize();

    /* - - - - - - - - - - - - A - - - - - - - - - - - - */
    for (int i = 0; i < nbNodesPerElementVelocity; i++){
      for (int j = 0; j < nbNodesPerElementVelocity; j++){
        /* P1 \oplus Bubble*/
        for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
          if (i < 3 and j < 3){
            /* P1 */  
            if (par.nbDofsPerNode()[0]*triangle[i]+comp >= m && par.nbDofsPerNode()[0]*triangle[i]+comp < n){
              double Mij = tria.mass(i,j);
              double Kij = tria.stiffness_symmetric(i,j);
              code = MatSetValue(M, par.nbDofsPerNode()[0]*triangle[i]+comp, par.nbDofsPerNode()[0]*triangle[j]+comp, 
                par.density()/par.timeStep() * Mij, ADD_VALUES); CHKERRQ(code);
              code = MatSetValue(B, par.nbDofsPerNode()[0]*triangle[i]+comp, par.nbDofsPerNode()[0]*triangle[j]+comp, 
                2.0 * par.viscosity() * Kij, ADD_VALUES); CHKERRQ(code);
            }

          } else if (i == 3 and j == 3) {
            /* Bubble */
            if (nbDofsVelP1 + par.nbDofsPerNode()[0]*triangle[i]+comp >= m 
              && nbDofsVelP1 + par.nbDofsPerNode()[0]*triangle[i]+comp < n){
              double Mij = tria.mass(i,j);
              double Kij = tria.stiffness_symmetric(i,j);
              code = MatSetValue(M, 
                nbDofsVelP1 + par.nbDofsPerNode()[0]*triangle[i]+comp, nbDofsVelP1 + par.nbDofsPerNode()[0]*triangle[j]+comp, 
                par.density()/par.timeStep() * Mij, ADD_VALUES); CHKERRQ(code);
              code = MatSetValue(B, 
                par.nbDofsPerNode()[0]*triangle[i]+comp, par.nbDofsPerNode()[0]*triangle[j]+comp, 
                2.0 * par.viscosity() * Kij, ADD_VALUES); CHKERRQ(code);
            }

          } else if (i < 3 and j == 3) {
            if (par.nbDofsPerNode()[0]*triangle[i]+comp >= m 
              && par.nbDofsPerNode()[0]*triangle[i]+comp < n){
              double Mij = tria.mass(i,j);
              double Kij = tria.stiffness_symmetric(i,j);
              code = MatSetValue(M, 
                par.nbDofsPerNode()[0]*triangle[i]+comp, nbDofsVelP1 + par.nbDofsPerNode()[0]*triangle[j]+comp, 
                par.density()/par.timeStep() * Mij, ADD_VALUES); CHKERRQ(code);
              code = MatSetValue(B, 
                par.nbDofsPerNode()[0]*triangle[i]+comp, par.nbDofsPerNode()[0]*triangle[j]+comp, 
                2.0 * par.viscosity() * Kij, ADD_VALUES); CHKERRQ(code);
            }

          } else if (i == 3 and j < 3) {
            /* Bubble */
            if (nbDofsVelP1 + par.nbDofsPerNode()[0]*triangle[i]+comp >= m 
              && nbDofsVelP1 + par.nbDofsPerNode()[0]*triangle[i]+comp < n){
              double Mij = tria.mass(i,j);
              double Kij = tria.stiffness_symmetric(i,j);
              code = MatSetValue(M, 
                nbDofsVelP1 + par.nbDofsPerNode()[0]*triangle[i]+comp, par.nbDofsPerNode()[0]*triangle[j]+comp, 
                par.density()/par.timeStep() * Mij, ADD_VALUES); CHKERRQ(code);
              code = MatSetValue(B, 
                par.nbDofsPerNode()[0]*triangle[i]+comp, par.nbDofsPerNode()[0]*triangle[j]+comp, 
                2.0 * par.viscosity() * Kij, ADD_VALUES); CHKERRQ(code);
            }
          }
        }
      }
    }

    /* - - - - - - - - - - - - (-B) - - - - - - - - - - - - */
    for (int i = 0; i < nbNodesPerElementVelocity; i++){
      for (int j = 0; j < nbNodesPerElementPressure; j++){
        for (int comp = 0; comp < nbDofsPerNode; comp++){
          if (i < 3 and j < 3){
            if (nbDofsPerNode*triangle[i]+comp >= m && nbDofsPerNode*triangle[i]+comp < n){
              double Bij = tria.mixed(i,j,comp);
              code = MatSetValue(B, 
                nbDofsPerNode*triangle[i]+comp, nbDofsVel + triangle[j], 
                -1.0*Bij, ADD_VALUES); CHKERRQ(code); /* -B*/
            }
          }
        }
      }
    }

    for (int i = 0; i < nbNodesPerElementVelocity; i++){
      for (int j = 0; j < nbNodesPerElementPressure; j++){
        for (int comp = 0; comp < nbDofsPerNode; comp++){
          if (i == 3 and j < 3){
            if (nbDofsVelP1 + nbDofsPerNode*triangle[i]+comp >= m && nbDofsPerNode*triangle[i]+comp < n){
              double Bij = tria.mixed(i,j,comp);
              code = MatSetValue(B, 
                nbDofsVelP1 + nbDofsPerNode*triangle[i]+comp, nbDofsVel + triangle[j], 
                -1.0*Bij, ADD_VALUES); CHKERRQ(code); /* -B*/
            }
          }
        }
      }
    }

    /* - - - - - - - - - - - - B^T - - - - - - - - - - - - */
    for (int i = 0; i < nbNodesPerElementPressure; i++){
      for (int j = 0; j < nbNodesPerElementVelocity; j++){
        for (int comp = 0; comp < nbDofsPerNode; comp++){
          if (i < 3 and j == 3){
            if (nbDofsVel + triangle[j] >= m && nbDofsVel + triangle[j] < n){
              /* div_phi_i phi_j */
              double Bij = tria.mixed(i,j,comp);
              code = MatSetValue(B, 
                nbDofsVel + triangle[j], nbDofsVelP1 + nbDofsPerNode * triangle[i] + comp, 
                1.0*Bij, ADD_VALUES); CHKERRQ(code); /* B^T */
            }
          }
        }
      }
    }

    for (int i = 0; i < nbNodesPerElementPressure; i++){
      for (int j = 0; j < nbNodesPerElementVelocity; j++){
        for (int comp = 0; comp < nbDofsPerNode; comp++){
          if (i < 3 and j < 3){
            if (nbDofsVel + triangle[j] >= m && nbDofsVel + triangle[j] < n){
              /* div_phi_i phi_j */
              double Bij = tria.mixed(i,j,comp);
              code = MatSetValue(B, nbDofsVel + triangle[j], nbDofsPerNode * triangle[i] + comp, 
                1.0*Bij, ADD_VALUES); CHKERRQ(code); /* B^T */
            }
          }
        }
      }
    }

    if (triaId % (geo.triangles()[0].size() / 5) == 0){
      if (world_rank == 0) cout << "  ++Elementary matrix for tria: " << triaId << "/" << geo.triangles()[0].size() - 1 << endl;
      if (world_rank == 0) cout << "  Size: " << tria.size() << endl;
      if (world_rank == 0) cout << "  Mass(1,1) = " << tria.mass(1,1) << endl;
      if (world_rank == 0) cout << "  Stiff(3,3) = " << tria.stiffness_symmetric(3,3) << endl;
      if (world_rank == 0) cout << "  Mixed(3,2,1) = " << tria.mixed(3,2,1) << endl;
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
  MatDiagonalSet(A, zeros(nbDofs), INSERT_VALUES);

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
    Vec rhs = zeros(nbDofs);
    code = MatMult(M, u0, rhs); CHKERRQ(code);
    if (world_rank == 0) cout << "Iteration: " << t/par.timeStep() << endl;

    /* Assembling convection */
    Mat C = mat(nbDofs, nbDofs);
//    Vec u0_seq = getSequential(u0);
//
//    if (world_rank == 0) cout << "\n - - - - - - - - - - - - - - - - - - - " << endl;
//    for (int triaId = 0; triaId < geo.triangles()[0].size(); triaId++){ 
//
//      vector<int> triangle = geo.triangles()[0][triaId];
//      /* get finite element coordinates */
//      vector<vector<double>> coordinates(nbNodesPerElement);
//      for (int i=0; i < nbNodesPerElement; i++){
//        coordinates[i] = geo.coordinates()[triangle[i]];
//      }
//      tria.setCoordinates(coordinates);
//
//    if (world_rank == 0) cout << "\n - - - - - - - - - - - - - - - - - - - " << endl;
//      /* loc2glob dof mapping */
//      vector<int> loc2glob(nbNodesPerElement*nbDofsPerNode);
//      for (int i = 0; i < nbNodesPerElement; i++){
//        for (int comp = 0; comp < nbDofsPerNode; comp++){
//          loc2glob[nbDofsPerNode*i+comp] = nbDofsPerNode*triangle[i]+comp;
//        }
//      }
//
//      /* Gather advection field on element */
//      vector<double> u_el(nbNodesPerElement*nbDofsPerNode);
//      code = VecGetValues(u0_seq, nbNodesPerElement*nbDofsPerNode, &loc2glob[0], &u_el[0]); CHKERR(code);
//
//      /* Set values from local to global */
//      for (int i = 0; i < nbNodesPerElement; i++){
//        for (int j = 0; j < nbNodesPerElement; j++){
//          vector<double> u_node_j(nbDofsPerNode);
//          for (int comp = 0; comp < nbDofsPerNode; comp++){
//            u_node_j[comp] = u_el[nbDofsPerNode*j+comp];
//          }
//          for (int comp = 0; comp < nbDofsPerNode; comp++){
//            if (nbDofsPerNode*triangle[i]+comp >= m && nbDofsPerNode*triangle[i]+comp < n){
//              code = MatSetValue(C, nbDofsPerNode*triangle[i]+comp, nbDofsPerNode*triangle[j]+comp, tria.advection(i, j, u_node_j), ADD_VALUES); CHKERRQ(code);
//            }
//          }
//        }
//      }
//      if (triaId % (geo.triangles()[0].size() / 5) == 0){
//        vector<double> u_node_j(nbDofsPerNode);
//        for (int comp = 0; comp < nbDofsPerNode; comp++){
//          u_node_j[comp] = u_el[nbDofsPerNode*0+comp];
//        }
//        if (world_rank == 0) cout << "  +++ Elementary matrix for triangle: " << triaId << "/" << geo.triangles()[0].size() - 1 << endl;
//        if (world_rank == 0) cout << "  |Jacobian|: " << tria.detJacobian() << endl;
//        if (world_rank == 0) cout << "  advection(0,1, u) = " << tria.advection(0,1,u_node_j) << endl;
//      }
//    }
    code = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
    code = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

    /* Add static part */
    MatAXPY(C, 1.0, A, DIFFERENT_NONZERO_PATTERN);

    boundary.block(rhs);
    boundary.block(C);

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
