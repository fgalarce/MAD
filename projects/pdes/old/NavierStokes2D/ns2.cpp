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
  vector<double> center(2, 0.0);
  double r = 1.0;
  double T = 0.5;
  double A = 60;
//  if (t < T/4){
//    bc[0] = A*(1.0 - x[1]*x[1]/(r*r) )*sin(2*PI/T*t);
//  } else {
    bc[0] = (1.0 / 0.41 / 0.41) * (1.2 * x[1] * (0.41 - x[1]));
//  }
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

  MasterTriangle2D tet; 
  tet.initialize(par, par.nbDofsPerNode()[0]);

  Boundary boundary;
  boundary.initialize(par, geo);

  int nbDofs = (par.nbDofsPerNode()[0] + par.nbDofsPerNode()[1])*io.nbVertices();
  int nbDofsPress = io.nbVertices();
  int nbDofsVel = par.nbDofsPerNode()[0]*io.nbVertices();
  int nbVertices = io.nbVertices();
  int nbDofsPerNode = par.nbDofsPerNode()[0];

  /* Save initial condition */
  Vec phi0 = zeros(nbDofs);
  Vec phi = zeros(nbDofs);
  io.writeState(phi0, 0.0);

  Mat M = mat(nbDofs, nbDofs); /* Mass matrix */
  Mat B = mat(nbDofs, nbDofs); /* Static matrix */

  int m,n;
  code = MatGetOwnershipRange(M, &m, &n); CHKERRQ(code);
  if (world_rank == 0) cout << "NS: Assembling discretization matrix." << endl;
  for (int tetraId = 0; tetraId < geo.triangles()[0].size(); tetraId++){ /* loop on tetra */

    vector<int> tetra = geo.triangles()[0][tetraId];
    /* get finite element coordinates */
    vector<vector<double>> coordinates(3);
    coordinates[0] = geo.coordinates()[tetra[0]];
    coordinates[1] = geo.coordinates()[tetra[1]];
    coordinates[2] = geo.coordinates()[tetra[2]];

    tet.setCoordinates(coordinates);
    tet.computeSize();

    if (tetraId % (geo.triangles()[0].size() / 5) == 0){
      if (world_rank == 0) cout << "  Elementary matrix for tetra: " << tetraId << "/" << geo.triangles()[0].size() - 1 << endl;
    }
    /* Assemble elementary matrices */
    for (int j = 0; j < 3; j++){
      for (int i = 0; i < 3; i++){
        for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
          if (par.nbDofsPerNode()[0]*tetra[i]+comp >= m && par.nbDofsPerNode()[0]*tetra[i]+comp < n){
            double Kij = tet.stiffness_symmetric(i,j);
//            double Kij = tet.stiffness(i,j);
            double Mij = par.density()/par.timeStep()*tet.mass(i,j);
            /* stiffness */
            code = MatSetValue(M, par.nbDofsPerNode()[0]*tetra[i]+comp, par.nbDofsPerNode()[0]*tetra[j]+comp, Mij, ADD_VALUES); CHKERRQ(code);
            MatSetValue(B, par.nbDofsPerNode()[0]*tetra[i]+comp, par.nbDofsPerNode()[0]*tetra[j]+comp, 
                  2.0 * par.viscosity() * Kij, ADD_VALUES); 
            /* div_phi_i phi_j */
            double Bij = tet.mixed(i,j,comp);
            code = MatSetValue(B, par.nbDofsPerNode()[0]*tetra[i]+comp, nbDofsVel + tetra[j], -1.0*Bij, ADD_VALUES); CHKERRQ(code);  /* -B*/
          }
          if (nbDofsVel + tetra[j] >= m && nbDofsVel + tetra[j] < n){
            double Bij = tet.mixed(i,j,comp);
            code = MatSetValue(B, nbDofsVel + tetra[j], par.nbDofsPerNode()[0]*tetra[i]+comp, 1.0*Bij, ADD_VALUES); CHKERRQ(code);  /* B^T */
          }
        }
        if (nbDofsVel + tetra[i] >= m && nbDofsVel + tetra[i] < n){
          double Kij = tet.stiffness(i,j);
          /* Brezzi-Pitkaranta stabilization for P1-P1 elements */
          MatSetValue(B, nbDofsVel + tetra[i], nbDofsVel + tetra[j], tet.size()*tet.size() * Kij, ADD_VALUES); 
        }
      }
    }
  }
  code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

//  saveMat(M, "M.txt");
//  saveMat(B, "B.txt");

  /* LHS static matrix */
  Mat A = mat(nbDofs, nbDofs);
double mmmmm;
MatNorm(M, NORM_FROBENIUS, &mmmmm);
if (world_rank == 0) cout << mmmmm << endl;
  code = MatDuplicate(M, MAT_COPY_VALUES, &A); CHKERRQ(code);

MatNorm(A, NORM_FROBENIUS, &mmmmm);
if (world_rank == 0) cout << mmmmm << endl;

MatNorm(M, NORM_FROBENIUS, &mmmmm);
if (world_rank == 0) cout << mmmmm << endl;

MatNorm(B, NORM_FROBENIUS, &mmmmm);
if (world_rank == 0) cout << mmmmm << endl;

  code = MatAXPY(A, 1.0, B, DIFFERENT_NONZERO_PATTERN); CHKERRQ(code);

MatNorm(A, NORM_FROBENIUS, &mmmmm);
if (world_rank == 0) cout << mmmmm << endl;

//  saveMat(A, "A.txt");
  /* Krylov solver */
  KSP ksp;
  code = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERR(code); 
  PetscBool reuse_prec = PETSC_TRUE;

  for (double t : range(0, par.nbIterations()*par.timeStep(), par.timeStep())){

    if (world_rank == 0) cout << "\n- - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
    if (world_rank == 0) cout << "NS: Solving Navier-Stokes equation for time: " << t << endl;

    /* block and assembly rhs */
    Vec rhs = vec(nbDofs);
    code = MatMult(M, phi0, rhs); CHKERRQ(code);

    /* Assemble time dependent matrix */
    Mat C = mat(nbDofs, nbDofs);
    code = MatGetOwnershipRange(C, &m, &n); CHKERRQ(code);

    Vec phi0_seq = getSequential(phi0);
    if (world_rank == 0) cout << "NS: norm phi0: " << norm(phi0_seq) << endl;

    for (int tetraId = 0; tetraId < geo.triangles()[0].size(); tetraId++){ /* loop on tetra */

      vector<int> tetra = geo.triangles()[0][tetraId];
      /* get finite element coordinates */
      vector<vector<double>> coordinates(3);
      coordinates[0] = geo.coordinates()[tetra[0]];
      coordinates[1] = geo.coordinates()[tetra[1]];
      coordinates[2] = geo.coordinates()[tetra[2]];

      tet.setCoordinates(coordinates);

      if (tetraId % (geo.triangles()[0].size() / 5) == 0){
        if (world_rank == 0) cout << "  Elementary matrix for tetra: " << tetraId << "/" << geo.triangles()[0].size() - 1 << endl;
      }

      /* loc2glob dof mapping */
      vector<int> loc2glob(6);
      for (int i = 0; i < 3; i++){
        for (int comp = 0; comp < 2; comp++){
          loc2glob[2*i+comp] = 2*tetra[i]+comp;
        }
      }

      /* Gather advection field on element */
      vector<double> u_el(6);
      code = VecGetValues(phi0, 6, &loc2glob[0], &u_el[0]); CHKERR(code);

      /* Set values from local to global */
      for (int j = 0; j < 3; j++){
        vector<double> u_node_j(2, 0.0);
        for (int comp = 0; comp < 2; comp++){
          u_node_j[comp] = u_el[2*j+comp];
        }
        for (int i = 0; i < 3; i++){
          for (int comp = 0; comp < 2; comp++){
            if (par.nbDofsPerNode()[0]*tetra[i]+comp >= m && par.nbDofsPerNode()[0]*tetra[i]+comp < n){
              double Aij = tet.advection(i, j, u_node_j);
//              code = MatSetValue(C, par.nbDofsPerNode()[0]*tetra[i]+comp, par.nbDofsPerNode()[0]*tetra[j]+comp, Aij, ADD_VALUES); CHKERRQ(code);
            }
          }
        }
      }
    }
    code = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
    code = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

    MatNorm(C, NORM_FROBENIUS, &mmmmm);
    if (world_rank == 0) cout << "norm C: " << mmmmm << endl;

    MatNorm(A, NORM_FROBENIUS, &mmmmm);
    if (world_rank == 0) cout << "norm A: " << mmmmm << endl;

    /* Add static part */
    MatAXPY(C, 1.0, A, DIFFERENT_NONZERO_PATTERN);

    MatNorm(C, NORM_FROBENIUS, &mmmmm);
    if (world_rank == 0) cout << "norm A + C: " << mmmmm << endl;

    /* block and assembly C */
    boundary.time(t);
    for (int i = 0; i < par.walls().size(); i++){
      boundary.Dirichlet(par.walls()[i], noslip, 0);
    }
    boundary.Dirichlet(par.inlet(), inlet_shape, 0);
    boundary.block(rhs);
    boundary.block(C);

    code = KSPSetOperators(ksp, C, C); CHKERR(code); 
    code = KSPSetType(ksp, KSPPREONLY); CHKERR(code); 
    PC pc;
    code = KSPGetPC(ksp, &pc); CHKERR(code); 
    code = PCSetType(pc, PCLU); CHKERR(code); 
    code = KSPSetReusePreconditioner(ksp, reuse_prec); CHKERR(code);
//    code = KSPSetTolerances(ksp, 1e-06, 1e-50, 20000, 20000);
    code = KSPSetFromOptions(ksp); CHKERR(code); 
//    int max_it;
//    double rtol, abstol, dtol;
//    KSPGetTolerances(ksp, &rtol, &abstol, &dtol, &max_it);
//    if (world_rank == 0) cout << max_it << endl;
//    if (world_rank == 0) cout << rtol << endl;
//    if (world_rank == 0) cout << abstol << endl;
//    if (world_rank == 0) cout << dtol << endl;

    if (world_rank == 0) cout << "KSP: Solving linear system" << endl;
    code = KSPSolve(ksp, rhs, phi); CHKERRQ(code); 
    int its;
    double rnorm;
    code = KSPGetResidualNorm(ksp,&rnorm); CHKERRQ(code);
    code = KSPGetIterationNumber(ksp,&its); CHKERRQ(code);
    if (world_rank == 0) cout << "KSP: " << its << " iterations. Residual norm = " << rnorm << endl;
  
    io.writeState(phi, t + par.timeStep());

    VecZeroEntries(phi0);
    code = VecAXPY(phi0, 1.0, phi); CHKERRQ(code);

    double norm_sol = norm(phi);
    if (world_rank == 0) cout << "KSP: Norm solution = " << norm_sol << endl;

    VecDestroy(&phi0_seq);
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
