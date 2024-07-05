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

#include <mad.hpp>

vector<double> inlet_shape(vector<double> x, double t, Parameters par){
  vector<double> bc(2, 0.0);
  double mean_velocity = par.inlet_u0();
  double shape = 6.0 / (0.41 * 0.41) * x[1] * ( 0.41 - x[1] );
  bc[0] = mean_velocity * shape;
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

  MasterElement fe; 
  fe.initialize(par, io.dimension());

  Boundary boundary;
  boundary.initialize(par, geo);

  int nbDofs = (par.nbDofsPerNode()[0] + par.nbDofsPerNode()[1])*io.nbVertices();
  int nbDofsPress = io.nbVertices();
  int nbDofsVel = par.nbDofsPerNode()[0]*io.nbVertices();
  int nbVertices = io.nbVertices();

  /* local values */
  int nbDofsPerNode = par.nbDofsPerNode()[0]; 
  int nbNodesPerElement = geo.dimension()+1; 

  /* initial condition */
  Vec phi0 = zeros(nbDofs);
  if (par.initial_condition_pre() != "none" && par.initial_condition_vel() != "none"){
    Vec vel0 = zeros(nbVertices*3);
    Vec pre0 = zeros(nbDofsPress);
    loadVec(vel0, par.initial_condition_vel());
    loadVec(pre0, par.initial_condition_pre());
    vector<double> stl_vec_3D = stl(vel0);
    vector<double> stl_vec(nbDofsVel);
    vector<double> stl_pre = stl(pre0);

    for (int i = 0; i < nbVertices; i++){
      for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){ /* this jumps the third coordinate always saved in ensight */
        stl_vec[par.nbDofsPerNode()[0]*i+comp] = stl_vec_3D[3*i+comp];
      }
    }
    
    int phi_low, phi_high;
    VecGetOwnershipRange(phi0, &phi_low, &phi_high);
    for (int i = 0; i < nbDofs; i++){
      if (i >= phi_low && i < phi_high){
        if (i < nbDofsVel){
          code = VecSetValue(phi0, i, stl_vec[i], INSERT_VALUES); CHKERRQ(code);
        } else {
          code = VecSetValue(phi0, i, stl_pre[i - nbDofsVel], INSERT_VALUES); CHKERRQ(code);
        }
      }
    }
    code = VecAssemblyBegin(phi0); CHKERRQ(code);
    code = VecAssemblyEnd(phi0); CHKERRQ(code);
  }
  io.writeState(phi0, 0.0);

  Vec phi = zeros(nbDofs);

  Mat M = mat(nbDofs, nbDofs); /* Mass matrix */
  Mat B = mat(nbDofs, nbDofs); /* Static matrix */

  int m,n;
  code = MatGetOwnershipRange(M, &m, &n); CHKERRQ(code);
  if (world_rank == 0) cout << "NS: Assembling discretization matrix." << endl;
  for (int feId = 0; feId < geo.elements.size(); feId++){ /* loop on elements */

    vector<int> simplex = geo.elements[feId];
    /* get finite element coordinates */
    vector<vector<double>> coordinates(nbNodesPerElement);
    for (int i = 0; i < nbNodesPerElement; i++){
      coordinates[i] = geo.coordinates()[simplex[i]];
    }

    fe.setCoordinates(coordinates);
    fe.computeSize();

    if (feId % (geo.elements.size() / 5) == 0){
      if (world_rank == 0) cout << "  Elementary matrix for tetra: " << feId << "/" << geo.elements.size() - 1 << endl;
    }
    /* Assemble elementary matrices */
    for (int j = 0; j < nbNodesPerElement; j++){
      for (int i = 0; i < nbNodesPerElement; i++){
        for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
          if (par.nbDofsPerNode()[0]*simplex[i]+comp >= m && par.nbDofsPerNode()[0]*simplex[i]+comp < n){
            double Kij = fe.stiffness(i,j);
            double Mij = par.density()/par.timeStep()*fe.mass(i,j);
            /* stiffness */
            MatSetValue(M, par.nbDofsPerNode()[0]*simplex[i]+comp, par.nbDofsPerNode()[0]*simplex[j]+comp, Mij, ADD_VALUES); 
            MatSetValue(B, par.nbDofsPerNode()[0]*simplex[i]+comp, par.nbDofsPerNode()[0]*simplex[j]+comp, par.viscosity() * Kij, ADD_VALUES); 
            /* div_phi_i phi_j */
            double Bij = fe.mixed(i,j,comp);
            MatSetValue(B, par.nbDofsPerNode()[0]*simplex[i]+comp, nbDofsVel + simplex[j], -1.0*Bij, ADD_VALUES);  /* -B*/
          }
          if (nbDofsVel + simplex[j] >= m && nbDofsVel + simplex[j] < n){
            double Bij = fe.mixed(i,j,comp);
            MatSetValue(B, nbDofsVel + simplex[j], par.nbDofsPerNode()[0]*simplex[i]+comp, 1.0*Bij, ADD_VALUES);  /* B^T */
          }
        }
        if (nbDofsVel + simplex[i] >= m && nbDofsVel + simplex[i] < n){
          double Kij = fe.stiffness(i,j);
          /* Brezzi-Pitkaranta stabilization for P1-P1 elements */
          MatSetValue(B, nbDofsVel + simplex[i], nbDofsVel + simplex[j], fe.size()*fe.size() * Kij, ADD_VALUES); 
        }
      }
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
  configureKSP(ksp, par.solver(), par.preconditioner(), par.monitorKSP(), par.use_solution_as_guess_KSP(), par.reuse_preconditioner());

  for (double t : range(0, par.nbIterations()*par.timeStep(), par.timeStep())){

    if (world_rank == 0) cout << "\n- - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
    if (world_rank == 0) cout << "NS: Solving Navier-Stokes equation for time: " << t << endl;

    /* block and assembly rhs */
    Vec rhs = vec(nbDofs);
    code = MatMult(M, phi0, rhs); CHKERRQ(code);

    /* Assemble time dependent matrix */
    Mat C = mat(nbDofs, nbDofs);

    Vec phi0_seq = getSequential(phi);

    if (world_rank == 0) cout << "NS: Assembling convection matrix." << endl;
    for (int feId = 0; feId < elements.size(); feId++){ /* loop on elements */

      vector<int> simplex = geo.elements[feId];
      /* get finite element coordinates */
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int i = 0; i < nbNodesPerElement; i++){
        coordinates[i] = geo.coordinates()[simplex[i]];
      }
      fe.setCoordinates(coordinates);
      fe.computeSize();

      if (feId % (geo.elements.size() / 5) == 0){
        if (world_rank == 0) cout << "  Elementary matrix for simplex: " << feId << "/" << geo.elements.size() - 1 << endl;
      }

      /* loc2glob dof mapping */
      vector<int> loc2glob(nbNodesPerElement*nbDofsPerNode);
      for (int i = 0; i < nbNodesPerElement; i++){
        for (int comp = 0; comp < nbDofsPerNode; comp++){
          loc2glob[nbDofsPerNode*i+comp] = nbDofsPerNode*simplex[i]+comp;
        }
      }

      /* Gather advection field on element */
      vector<double> u_el(nbNodesPerElement*nbDofsPerNode);
      code = VecGetValues(phi0_seq, nbNodesPerElement*nbDofsPerNode, &loc2glob[0], &u_el[0]); CHKERR(code);

      /* Set values from local to global */
      for (int j = 0; j < nbNodesPerElement; j++){
        vector<double> u_node_j(par.nbDofsPerNode()[0]);
        for (int i = 0; i < nbNodesPerElement; i++){
          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
            u_node_j[comp] = u_el[par.nbDofsPerNode()[0]*j+comp];
          }
          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
            if (par.nbDofsPerNode()[0]*simplex[i]+comp >= m && par.nbDofsPerNode()[0]*simplex[i]+comp < n){
              double Aij = par.density()*fe.advection(i, j, u_node_j);
              double SUPG = fe.size()*fe.size() * fe.supg(i, j, u_node_j);
              MatSetValue(C, par.nbDofsPerNode()[0]*simplex[i]+comp, par.nbDofsPerNode()[0]*simplex[j]+comp, Aij + SUPG, ADD_VALUES);
            }
          }
        }
      }
    }
    code = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
    code = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

    /* Add static part */
    MatAXPY(C, 1.0, A, DIFFERENT_NONZERO_PATTERN);

    /* block and assembly C */
    boundary.time(t);
    boundary.Dirichlet(par.inlet(), inlet_shape, 0);
    for (int i = 0; i < par.walls().size(); i++){
      boundary.Dirichlet(par.walls()[i], noslip, 0);
    }
    boundary.block(rhs);
    boundary.block(C);

    code = KSPSetOperators(ksp, C, C); CHKERR(code); 
    code = KSPSetUp(ksp); CHKERRQ(code);

    if (world_rank == 0) cout << "KSP: Solving linear system" << endl;
    code = KSPSolve(ksp, rhs, phi); CHKERRQ(code); 

    if (!par.monitorKSP()){
      int its;
      double rnorm;
      code = KSPGetResidualNorm(ksp,&rnorm); CHKERRQ(code);
      code = KSPGetIterationNumber(ksp,&its); CHKERRQ(code);
      if (world_rank == 0) cout << "KSP: " << its << " iterations. Residual norm = " << rnorm << endl;
    }

    if (par.backflowStab()){
      boundary.backflow(phi,0);
    }

    io.writeState(phi, t + par.timeStep());

    VecZeroEntries(phi0);
    code = VecAXPY(phi0, 1.0, phi); CHKERRQ(code);

    double norm_sol = norm(phi);
    if (world_rank == 0) cout << "KSP: Norm solution = " << norm_sol << endl;

//    VecDestroy(&phi0_seq);
    VecDestroy(&rhs);
//    MatDestroy(&C);
  }

  KSPDestroy(&ksp);
  MatDestroy(&A);
  VecDestroy(&phi);
  VecDestroy(&phi0);
  MatDestroy(&B);
  MatDestroy(&M);
  SlepcFinalize();
}
