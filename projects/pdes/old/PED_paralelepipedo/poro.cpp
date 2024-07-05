/*=============================================================================
  This file is part of the code MAD (or MAD)
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
#include "bcs.hpp"

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

  Calculus cal;
  cal.initialize(par, geo, io);

  MasterElement fe; 
  fe.initialize(par, io.dimension());

  Boundary boundary;
  boundary.initialize(par, geo);

  int nbDofsPerNode = par.nbDofsPerNode()[0]; 
  int nbVertices = io.nbVertices();
  int nbDofsSol = nbVertices * nbDofsPerNode;
  int nbDofsFlu = nbVertices * nbDofsPerNode;
  int nbDofsPress = nbVertices;
  int nbDofs = nbDofsSol + nbDofsFlu + nbDofsPress;

  /* local values */
  int nbNodesPerElement = geo.dimension()+1; 

  /* initial condition */
  Vec u00 = zeros(nbDofs);
  Vec u0 = zeros(nbDofs);
  Vec u = zeros(nbDofs);

//  io.writeState(u0, 0.0);

  Mat M_sol_1 = mat(nbDofs, nbDofs); /* Mass matrix */
  Mat M_sol_2 = mat(nbDofs, nbDofs); /* Mass matrix */
  Mat A = mat(nbDofs, nbDofs); /* Static matrix */

  double K_darcy =  par.permeability()/par.viscosity();
  if (world_rank == 0) cout << "BIOT: permeability/viscosity = " << K_darcy << endl;

  int m,n;
  code = MatGetOwnershipRange(A, &m, &n); CHKERRQ(code);
  if (world_rank == 0) cout << "PORO: Assembling discretization matrix." << endl;
  for (int partId = 0; partId < geo.elements.size(); partId++){
    for (int feId = 0; feId < geo.elements[partId].size(); feId++){ /* loop on elements */

      vector<int> simplex = geo.elements[partId][feId];
      /* get finite element coordinates */
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int i = 0; i < nbNodesPerElement; i++){
        coordinates[i] = geo.coordinates()[simplex[i]];
      }

      fe.setCoordinates(coordinates);
      fe.computeSize();

      if (feId % (geo.elements[partId].size() / par.verbose()) == 0){
        if (world_rank == 0) cout << "  Elementary matrix for tetra: " << feId << "/" << geo.elements[partId].size() - 1 << endl;
      }
      /* Assemble elementary matrices */
      for (int j = 0; j < nbNodesPerElement; j++){
          
        /* Compute Lame coefficients */
        double mu, lambda;
        if (      coordinates[j][0] >= par.x_range()[0] && coordinates[j][0] <= par.x_range()[1]
               && coordinates[j][1] >= par.y_range()[0] && coordinates[j][1] <= par.y_range()[1]
               && coordinates[j][2] >= par.z_range()[0] && coordinates[j][2] <= par.z_range()[1] ){
            mu=60000; 
            lambda=45000;

//          lambda = par.lambda()[0];
//          mu = par.mu()[0];
//          lambda = par.youngModulus()[0]/(2.0*(1.0+par.poissonRatio()[0]));
//          mu = par.youngModulus()[0] * par.poissonRatio()[0] / ((1.0 + par.poissonRatio()[0])*(1.0 - 2.0*par.poissonRatio()[0]));
        } else {
            mu=30000; 
            lambda=90000;

//          lambda = par.lambda()[1];
//          mu = par.mu()[1];
//          lambda = par.youngModulus()[1]/(2.0*(1.0+par.poissonRatio()[1]));
//          mu = par.youngModulus()[1] * par.poissonRatio()[1] / ((1.0 + par.poissonRatio()[1])*(1.0 - 2.0*par.poissonRatio()[1]));
        }

        for (int i = 0; i < nbNodesPerElement; i++){
          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){

            /* Newton law for solid  */
            if (nbDofsPerNode*simplex[i]+comp >= m && nbDofsPerNode*simplex[i]+comp < n){
              double Kij = fe.stiffness_symmetric(i,j);
              double div_div_ij = fe.div_div(i,j,comp);
              double Mij = par.density() * fe.mass(i,j) / (par.timeStep()*par.timeStep());
              code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp, mu*Kij + lambda*div_div_ij, ADD_VALUES); CHKERRQ(code);
              code = MatSetValue(M_sol_1, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp, Mij, ADD_VALUES); CHKERRQ(code);
              /* div_phi_i phi_j */
              double Bij = fe.mixed(i,j,comp);
              code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp, nbDofsSol + nbDofsFlu + simplex[j], -1.0*Bij, ADD_VALUES); CHKERRQ(code);  /* -B*/ 
            }

            /* Darcy law for fluid */
            if (nbDofsSol + nbDofsPerNode*simplex[i]+comp >= m && nbDofsSol + nbDofsPerNode*simplex[i]+comp < n){
              double Mij = fe.mass(i,j);
              code = MatSetValue(A, nbDofsSol + nbDofsPerNode*simplex[i]+comp, nbDofsSol + nbDofsPerNode*simplex[j]+comp, K_darcy * Mij, ADD_VALUES);
              double Bij = fe.mixed(i,j,comp);
              code = MatSetValue(A, nbDofsSol + nbDofsPerNode*simplex[i]+comp, nbDofsSol + nbDofsFlu + simplex[j], -1.0*Bij, ADD_VALUES);
            }

            /* Mass conservation */
            if (nbDofsSol + nbDofsFlu + simplex[j] >= m && nbDofsSol + nbDofsFlu + simplex[j] < n){
              double Bij = fe.mixed(i,j,comp);
              MatSetValue(M_sol_2, nbDofsSol + nbDofsFlu + simplex[j], nbDofsPerNode*simplex[i]+comp, -1.0*Bij/par.timeStep(), ADD_VALUES);  /* B^T */
              MatSetValue(A, nbDofsSol + nbDofsFlu + simplex[j], nbDofsSol + nbDofsPerNode*simplex[i]+comp, -1.0*Bij, ADD_VALUES);  /* B^T */
              double Kij = fe.stiffness(i,j);
              /* Brezzi-Pitkaranta stabilization for P1-P1 geo.elements[partId] */
              code = MatSetValue(A, nbDofsSol + nbDofsFlu + simplex[i], nbDofsSol + nbDofsFlu + simplex[j], fe.size()*fe.size() * Kij, ADD_VALUES); CHKERRQ(code);
            }
          }
        }
      }
    }
  }
  code = MatAssemblyBegin(M_sol_1, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(M_sol_1, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyBegin(M_sol_2, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(M_sol_2, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

  /* LHS static matrix */
  Mat C = mat(nbDofs, nbDofs);
  code = MatDuplicate(A, MAT_COPY_VALUES, &C); CHKERRQ(code);
  code = MatAXPY(C, 1.0, M_sol_1, DIFFERENT_NONZERO_PATTERN); CHKERRQ(code);
  code = MatAXPY(C, 1.0, M_sol_2, DIFFERENT_NONZERO_PATTERN); CHKERRQ(code);

  /* Boundary conditions */
  boundary.Dirichlet(par.inlet(), inlet, 0);
  boundary.Dirichlet(par.outlet(), outlet, 0);
  boundary.Dirichlet(par.outlet(), outlet, 1);

  /* Krylov solver */
  KSP ksp;
  configureKSP(ksp, par.solver(), par.preconditioner(), par.monitorKSP(), par.use_solution_as_guess_KSP(), par.reuse_preconditioner());
  boundary.block(C);

  if (world_rank == 0) cout << "Elasticity: Seting operators" << endl;
  code = KSPSetOperators(ksp, C, C); CHKERR(code); 
  code = KSPSetUp(ksp); CHKERRQ(code);

  ofstream file_norm_solution(par.dirResults() + "/norm_solution.txt");
  
  for (double t : range(0, par.nbIterations()*par.timeStep(), par.timeStep())){

    if (world_rank == 0) cout << "\n- - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
    if (world_rank == 0) cout << "Poro: Solving poro-elasticity equation for time: " << t << endl;

    /* Compute velocity to couple solid and fluid */
    Vec backward_velocity = zeros(nbDofs); 
    code = VecAXPY(backward_velocity,  1.0, u0); CHKERRQ(code);
    code = VecAXPY(backward_velocity, -1.0, u00); CHKERRQ(code);
    code = VecScale(backward_velocity, 1.0/par.timeStep()); CHKERRQ(code);

    /* RHS fluid */
    boundary.time(t); 
    Vec b_sol2 = zeros(nbDofs);
    code = MatMult(M_sol_2, u0, b_sol2); CHKERRQ(code);
    double Mu_sol2_norm = norm(b_sol2);
    if (world_rank == 0) cout << "Elasticity: Norm b_solid_2 = " << Mu_sol2_norm << endl;

    /* RHS solid */
    code = VecScale(u0, 2.0); CHKERRQ(code);
    code = VecAXPY(u0, -1.0, u00); CHKERRQ(code);
    Vec b = zeros(nbDofs);
    code = MatMult(M_sol_1, u0, b); CHKERRQ(code);
    double Mu_sol_norm = norm(b);
    if (world_rank == 0) cout << "Elasticity: Norm b_solid_1 = " << Mu_sol_norm << endl;

    code = VecAXPY(b, 1.0, b_sol2); CHKERRQ(code);

    /* Boundary conditions */
    boundary.Dirichlet(par.inlet(), inlet, 0);
    boundary.Dirichlet(par.outlet(), outlet, 0);
    boundary.Dirichlet(par.outlet(), outlet, 1);
   
    /* block and set up system of equations */ 
    boundary.block(b);
    if (world_rank == 0) cout << "Elasticity: Solving linear system" << endl;
    code = KSPSolve(ksp, b, u); CHKERRQ(code); 
    int its;
    double rnorm;
    code = KSPGetResidualNorm(ksp, &rnorm); CHKERRQ(code);
    code = KSPGetIterationNumber(ksp, &its); CHKERRQ(code);

    if (world_rank == 0) cout << "Elasticity: " << its << " iterations. Residual norm = " << rnorm << endl;
    
    io.writeState(u, t + par.timeStep());
       
    /* Update solutions for centered finite differences */ 
    if (t > 0){
      VecZeroEntries(u00); 
      VecAXPY(u00, 1.0, u0);
    }
    code = VecZeroEntries(u0); CHKERRQ(code);
    code = VecAXPY(u0, 1.0, u); CHKERRQ(code);

    double norm_sol = norm(u);
    if (world_rank == 0) cout << "Elasticity: Norm solution = " << norm_sol << endl;
    file_norm_solution << norm_sol << endl; 

  }
  file_norm_solution.close();

  KSPDestroy(&ksp);
  MatDestroy(&A);
  VecDestroy(&u);
  VecDestroy(&u0);
  VecDestroy(&u00);
  MatDestroy(&C);
  MatDestroy(&M_sol_1);
  MatDestroy(&M_sol_2);
  par.finalize();
  SlepcFinalize();
}
