/*=============================================================================
  This file is part of the code MAD
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2021,
    
     Felipe Galarce

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

/* function pointers are passed as argument to BC class */
vector<double> inlet(vector<double> x, double t){
  vector<double> bc(2, 0.0);
  bc[0] = 1.0;
  return bc;
}

vector<double> outlet(vector<double> x, double t){
  vector<double> bc(2, 0.0);
  return bc;
}

vector<double> wall(vector<double> x, double t){
  vector<double> bc(1, 0.0);
  return bc;
}

void matSet(Mat A, double i, double j, double Aij){
  int m,n;
  MatGetOwnershipRange(A, &m, &n); 
  if (i >= m && i < n){
    MatSetValue(A, i, j, Aij, ADD_VALUES);
  }
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

  Calculus cal;
  cal.initialize(par, geo, io);

  MasterTriangle2D tria; 
  int nbDofsPerNode = par.nbDofsPerNode()[0]; // scalar field
  tria.initialize(par, nbDofsPerNode);

  int nbDofs_uSol = io.nbVertices() * nbDofsPerNode;
  int nbDofs_uFlu = io.nbVertices() * nbDofsPerNode;
  int nbDofsPress = io.nbVertices();
  int nbDofs = nbDofs_uSol + nbDofs_uFlu + nbDofsPress;

  Mat A = mat(nbDofs, nbDofs); /* Static matrix */
  Mat M = mat(nbDofs, nbDofs); /* Mass matrix */


  Vec u0 = zeros(nbDofs);
  Vec u = zeros(nbDofs);

  /* Boundary conditions */
  BoundaryConditions bc;
  bc.initialize(par, geo);
  bc.Dirichlet(par.inlet(), inlet, 0);
  bc.Dirichlet(par.inlet(), outlet, 1);
  bc.Dirichlet(par.outlet(), outlet, 0);
  bc.Dirichlet(par.outlet(), outlet, 1);
  bc.block(u0);

  io.writeState(u0, 0.0);

//  /* Assemble mu */

  if (world_rank == 0) cout << "BIOT: Assembling discretization matrix" << endl;

  for (int triaId = 0; triaId < geo.triangles()[0].size(); triaId++){ /* loop on tria */
    vector<int> triangle = geo.triangles()[0][triaId];

    /* get finite element coordinates */
    vector<vector<double>> coordinates(3);
    coordinates[0] = geo.coordinates()[triangle[0]];
    coordinates[1] = geo.coordinates()[triangle[1]];
    coordinates[2] = geo.coordinates()[triangle[2]];

    tria.setCoordinates(coordinates);
    tria.computeSize();
  
    if (triaId % (geo.triangles()[0].size() / 5) == 0){
      if (world_rank == 0) cout << "  \n* * * * * * * * * * * * * * * " << endl;
      if (world_rank == 0) cout << "  Elementary matrix for tria: " << triaId << "/" << geo.triangles()[0].size() - 1 << endl;
      if (world_rank == 0) cout << "  Coordinates: " << endl << coordinates;
      if (world_rank == 0) cout << "  Labels: " << triangle << endl; 
      if (world_rank == 0) cout << "  Size: " << tria.size() << endl;
      if (world_rank == 0) cout << "  |Jacobian| = " << tria.detJacobian() << endl;
    }

    /* Assemble elementary matrices */
    for (int i = 0; i < 3; i++){
      for (int j = 0; j < 3; j++){

        double Mij = tria.mass(i,j);
        double Kij = tria.stiffness(i,j);

        for (int comp = 0; comp < nbDofsPerNode; comp++){
          double Bij = tria.mixed(i,j,comp);
          matSet(A, nbDofsPerNode*triangle[i]+comp, nbDofsPerNode*triangle[j]+comp, par.viscosity() * Kij);
          matSet(A, nbDofs_uSol + nbDofsPerNode*triangle[i]+comp, nbDofs_uSol + nbDofsPerNode*triangle[j]+comp, 1.0/par.porosity() * Mij);
          matSet(A, nbDofsPerNode * triangle[i] + comp, nbDofs_uSol + nbDofs_uFlu + triangle[j], -1.0*Bij);
          matSet(A, nbDofs_uSol + nbDofsPerNode * triangle[i] + comp, nbDofs_uSol + nbDofs_uFlu + triangle[j], -1.0*Bij);
          matSet(A, nbDofs_uSol + nbDofs_uFlu + triangle[j], nbDofsPerNode * triangle[i] + comp, par.timeStep() * Bij);
          matSet(A, nbDofs_uSol + nbDofs_uFlu + triangle[j], nbDofs_uSol + nbDofsPerNode * triangle[i] + comp, Bij);
        }
        matSet(M, nbDofs_uSol + nbDofs_uFlu + triangle[i], triangle[j], Mij);
        matSet(A, nbDofs_uSol + nbDofs_uFlu + triangle[i], nbDofs_uSol + nbDofs_uFlu + triangle[j], tria.size()*tria.size() * Kij);
      }
    }
  }
  code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

  bc.block(A);

  /* Krylov solver */
  KSP ksp;
  code = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERR(code); 
  PetscBool reuse_prec = PETSC_TRUE;
  code = KSPSetOperators(ksp, A, A); CHKERR(code); 
  code = KSPSetType(ksp, KSPGMRES); CHKERR(code); 
  PC pc;
  code = KSPGetPC(ksp, &pc); CHKERR(code); 
  code = PCSetType(pc, PCASM); CHKERR(code); 
  code = KSPSetReusePreconditioner(ksp, reuse_prec); CHKERR(code);
  code = KSPSetFromOptions(ksp); CHKERR(code); 

  ofstream ksp_output(par.dirResults() + "/ksp.log");

  for (double t : range(0, par.nbIterations()*par.timeStep(), par.timeStep())){

    if (world_rank == 0) cout << "\n - - - - - - - - - - - - - - - - - - - " << endl;
    if (world_rank == 0) cout << "Solving Poro-elasticity equations for time: " << t << endl;
    if (world_rank == 0) cout << "Iteration: " << t/par.timeStep() << endl;

    /* compute div(u0) */
    Vec uSol_vec = vec(nbDofs_uSol);
    double * uSolArr = new double[nbDofs_uSol];
    code = VecGetValues(u0, nbDofs_uSol, &range(nbDofs_uSol)[0], uSolArr); CHKERRQ(code);
    code = VecSetValues(uSol_vec, nbDofs_uSol, &range(nbDofs_uSol)[0], uSolArr, INSERT_VALUES); CHKERRQ(code);
    code = VecAssemblyBegin(uSol_vec); CHKERRQ(code);
    code = VecAssemblyEnd(uSol_vec); CHKERRQ(code);

    Vec div_u0 = cal.divergence(uSol_vec);
    double * divU0arr = new double[nbDofsPress];
    Vec b0 = zeros(nbDofs);
    code = VecGetValues(div_u0, nbDofsPress, &range(nbDofsPress)[0], divU0arr); CHKERRQ(code);
    code = VecSetValues(b0, nbDofsPress, &range(nbDofs_uSol+nbDofs_uFlu, nbDofs_uSol+nbDofs_uFlu+nbDofsPress)[0], divU0arr, INSERT_VALUES); CHKERRQ(code);
    code = VecAssemblyBegin(b0); CHKERRQ(code);
    code = VecAssemblyEnd(b0); CHKERRQ(code);
    cout << "BIOT: norm(div u0) = " << norm(div_u0) << endl;

    /* block and assembly b */
    Vec b = zeros(nbDofs);
    code = MatMult(M, b0, b); CHKERRQ(code);
    bc.block(b);

    if (world_rank == 0) cout << "BIOT: Solving linear system" << endl;
    code = KSPSolve(ksp, b, u); CHKERRQ(code); 
    int its;
    double rnorm;
    code = KSPGetResidualNorm(ksp, &rnorm); CHKERRQ(code);
    code = KSPGetIterationNumber(ksp, &its); CHKERRQ(code);
    if (world_rank == 0) cout << "KSP: " << its << " iterations. Residual norm = " << rnorm << endl;
    ksp_output << its << " " << rnorm << endl;

    io.writeState(u, t);

    VecZeroEntries(u0);
    code = VecAXPY(u0, 1.0, u); CHKERRQ(code);

    double norm_sol = norm(u);
    if (world_rank == 0) cout << "  Norm solution = " << norm_sol << endl;

    VecDestroy(&b);
  }
  ksp_output.close();

  KSPDestroy(&ksp);
  MatDestroy(&A);
  VecDestroy(&u);
  VecDestroy(&u0);
  MatDestroy(&M);
  par.finalize();
  SlepcFinalize();
}
