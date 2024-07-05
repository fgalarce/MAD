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
vector<double> inlet(vector<double> x, double t, Parameters par){
  vector<double> bc(2, 0.0);
  bc[0] = 1.0;
  return bc;
}

vector<double> outlet(vector<double> x, double t, Parameters par){
  vector<double> bc(2, 0.0);
  return bc;
}

vector<double> wall(vector<double> x, double t, Parameters par){
  vector<double> bc(1, 0.0);
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

  Calculus cal;
  cal.initialize(par, geo, io);

  MasterElement fe; 
  fe.initialize(par, io.dimension());

  /* local values */
  int nbDofsPerNode = par.nbDofsPerNode()[0]; 
  int nbNodesPerElement = geo.dimension()+1; 

  /* global values */
  int nbDofs_uSol = io.nbVertices() * nbDofsPerNode;
  int nbDofs_uFlu = io.nbVertices() * nbDofsPerNode;
  int nbDofsPress = io.nbVertices();
  int nbDofs = nbDofs_uSol + nbDofs_uFlu + nbDofsPress;
  int nbVertices = io.nbVertices();

  Mat A = mat(nbDofs, nbDofs); /* Static matrix */
  Mat M = mat(nbDofs, nbDofs); /* Mass matrix */

  Vec u0 = zeros(nbDofs);
  Vec u = zeros(nbDofs);

  /* Boundaries handler */
  Boundary boundary;
  boundary.initialize(par, geo);
  boundary.Dirichlet(par.inlet(), inlet, 0);
  boundary.Dirichlet(par.inlet(), outlet, 1);
  boundary.Dirichlet(par.outlet(), outlet, 0);
  boundary.Dirichlet(par.outlet(), outlet, 1);
  boundary.block(u0);

  io.writeState(u0, 0.0);

  vector<vector<int>> elements;
  if (geo.dimension() == 2){
    elements = geo.triangles()[0];
  } else {
    elements = geo.tetrahedron()[0];
  }

  int m,n;
  code = MatGetOwnershipRange(M, &m, &n); CHKERRQ(code);

  /* Compute Lame coefficients */
  double lambda = par.youngModulus()/(2.0*(1.0+par.poissonRatio()));
  double mu = par.youngModulus() * par.poissonRatio() / ((1.0 + par.poissonRatio())*(1.0 - 2.0*par.poissonRatio()));
  if (world_rank == 0) cout << "Elasticity: Lame parameters. mu = " << mu << ". lambda = " << lambda << ". " << endl;

  if (world_rank == 0) cout << "BIOT: Assembling discretization matrix" << endl;

  for (int feId = 0; feId < elements.size(); feId++){ /* loop on tria */

    vector<int> simplex = elements[feId];
    /* get finite element coordinates */
    vector<vector<double>> coordinates(nbNodesPerElement);
    for (int i = 0; i < nbNodesPerElement; i++){
      coordinates[i] = geo.coordinates()[simplex[i]];
    }

    fe.setCoordinates(coordinates);
    fe.computeSize();
  
    if (feId % (elements.size() / 5) == 0){
      if (world_rank == 0) cout << "  --> Elementary matrix for fe: " << feId << "/" << elements.size() - 1 << endl;
      if (world_rank == 0) cout << "  Coordinates: " << endl << coordinates;
      if (world_rank == 0) cout << "  Size: " << fe.size() << endl;
      if (world_rank == 0) cout << "  |Jacobian| = " << fe.detJacobian() << endl;
    }

    /* Assemble elementary matrices */
    for (int i = 0; i < nbNodesPerElement; i++){
      for (int j = 0; j < nbNodesPerElement; j++){

        double Mij = fe.mass(i,j);
        double Kij_BP = fe.stiffness(i,j);

        for (int comp = 0; comp < nbDofsPerNode; comp++){
          if (nbDofsPerNode*simplex[i]+comp >= m && nbDofsPerNode*simplex[i]+comp < m){
            double Kij = fe.stiffness_symmetric(i,j);
            code = MatSetValue(A, 
                  nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp, 
                  par.mu() * Kij); CHKERRQ(code);
            code = MatSetValue(A, 
                  nbDofsPerNode * simplex[i] + comp, nbDofs_uSol + nbDofs_uFlu + simplex[j], 
                  -1.0*Bij);

          }
          double Bij = fe.mixed(i,j,comp);
          code = MatSetValue(A, 
                nbDofs_uSol + nbDofsPerNode*simplex[i]+comp, nbDofs_uSol + nbDofsPerNode*simplex[j]+comp, 
                1.0/par.porosity() * Mij);
          code = MatSetValue(A, 
                nbDofs_uSol + nbDofsPerNode * simplex[i] + comp, nbDofs_uSol + nbDofs_uFlu + simplex[j], 
                -1.0*Bij);
          code = MatSetValue(A, 
                nbDofs_uSol + nbDofs_uFlu + simplex[j], nbDofsPerNode * simplex[i] + comp, 
                par.timeStep() * Bij);
          code = MatSetValue(A, 
                nbDofs_uSol + nbDofs_uFlu + simplex[j], nbDofs_uSol + nbDofsPerNode * simplex[i] + comp, 
                Bij);
        }
        code = MatSetValue(M, 
                nbDofs_uSol + nbDofs_uFlu + simplex[i], simplex[j], 
                Mij);
        code = MatSetValue(A, 
                nbDofs_uSol + nbDofs_uFlu + simplex[i], nbDofs_uSol + nbDofs_uFlu + simplex[j], 
                fe.size()*fe.size() * Kij_BP);
      }
    }
  }
  code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

  boundary.block(A);

  /* Krylov solver */
  KSP ksp;
  configureKSP(ksp, par.solver(), par.preconditioner(), par.monitorKSP(), par.use_solution_as_guess_KSP(), par.reuse_preconditioner());
  code = KSPSetOperators(ksp, A, A); CHKERR(code); 
  code = KSPSetUp(ksp); CHKERRQ(code);

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
    boundary.block(b);

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
