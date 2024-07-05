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

vector<double> zero(vector<double> x, double t, Parameters par){
  vector<double> bc(1, 0.0);
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

  int nbDofsPerNode = par.nbDofsPerNode()[0]; 
  int nbNodesPerElement = io.dimension()+1; 

  int nbDofs_uSol = io.nbVertices() * nbDofsPerNode;
  int nbDofs_uFlu = io.nbVertices() * nbDofsPerNode;
  int nbDofsPress = io.nbVertices();
  int nbDofs = nbDofs_uSol + nbDofs_uFlu + nbDofsPress;

  Mat A = mat(nbDofs, nbDofs); /* Static matrix */
  Mat M = mat(nbDofs, nbDofs); /* Mass matrix */

  Vec u0 = zeros(nbDofs);
  Vec u = zeros(nbDofs);

  /* Boundary conditions */
  Boundary boundary;
  boundary.initialize(par, geo);
  boundary.Dirichlet(par.inlet(), inlet, 0);
  boundary.Dirichlet(par.inlet(), inlet, 1);
  boundary.Dirichlet(par.outlet(), outlet, 0);
  boundary.Dirichlet(par.outlet(), outlet, 1);
  for (int i : par.walls()){
    boundary.DirichletComp(0, i, wall, 0);
    boundary.DirichletComp(1, i, wall, 0);
    boundary.DirichletComp(0, i, wall, 1);
    boundary.DirichletComp(1, i, wall, 1);
  }
  boundary.Dirichlet(par.outlet(), zero, 2);

  io.writeState(u0, 0.0);

  vector<vector<int>> elements;
  if (geo.dimension() == 2){
    elements = geo.triangles()[0];
  } else {
    elements = geo.tetrahedron()[0];
  }

  int m,n;
  code = MatGetOwnershipRange(A, &m, &n); CHKERRQ(code);

  /* Compute Lame coefficients */
  double lambda = par.youngModulus()/(2.0*(1.0+par.poissonRatio()));
  double mu = par.youngModulus() * par.poissonRatio() / ((1.0 + par.poissonRatio())*(1.0 - 2.0*par.poissonRatio()));
  if (world_rank == 0) cout << "BIOT: Lame parameters. mu = " << mu << ". lambda = " << lambda << ". " << endl;

  double K_darcy =  par.permeability()*par.porosity()/par.viscosity();
  if (world_rank == 0) cout << "BIOT: permeability*porosity/viscosity = " << K_darcy << endl;

  if (world_rank == 0) cout << "BIOT: Assembling discretization matrix" << endl;

  for (int feId = 0; feId < elements.size(); feId++){ /* loop on fe */
    vector<int> elementNodes = elements[feId];

    /* get finite element coordinates */
    vector<vector<double>> coordinates(nbNodesPerElement);
    for (int nodeId = 0; nodeId < nbNodesPerElement; nodeId++){
      coordinates[nodeId] = geo.coordinates()[elementNodes[nodeId]];
    }
    
    fe.setCoordinates(coordinates);
    fe.computeSize();
  
    if (feId % (elements.size() / 5) == 0){
      if (world_rank == 0) cout << "  Elementary matrix for tetra: " << feId << "/" << elements.size() - 1 << endl;
      if (world_rank == 0) cout << "  Size: " << fe.size() << endl;
      if (world_rank == 0) cout << "  detJac: " << fe.detJacobian() << endl;
    }
    /* Assemble elementary matrices */
    for (int i = 0; i < nbNodesPerElement; i++){
      for (int j = 0; j < nbNodesPerElement; j++){
        for (int comp = 0; comp < nbDofsPerNode; comp++){
          if (nbDofsPerNode*elementNodes[i]+comp >= m && nbDofsPerNode*elementNodes[i]+comp < n){
            double Bij = fe.mixed(i,j,comp);
            double Kij = fe.stiffness_symmetric(i,j);
            double div_div_ij = fe.div_div(i,j,comp);
            MatSetValue(A, nbDofsPerNode*elementNodes[i]+comp, nbDofsPerNode*elementNodes[j]+comp, mu*Kij + lambda*div_div_ij, ADD_VALUES); 
            MatSetValue(A, nbDofsPerNode*elementNodes[i]+comp, nbDofs_uSol + nbDofs_uFlu + elementNodes[j], -1.0*Bij, ADD_VALUES);
          }
          if (nbDofs_uSol + nbDofsPerNode*elementNodes[i]+comp >= m && nbDofs_uSol + nbDofsPerNode*elementNodes[i]+comp < n){
            double Bij = fe.mixed(i,j,comp);
            double Mij = fe.mass(i,j);
            MatSetValue(A, nbDofs_uSol + nbDofsPerNode*elementNodes[i]+comp, nbDofs_uSol + nbDofsPerNode*elementNodes[j]+comp,  K_darcy * Mij, ADD_VALUES);
            MatSetValue(A, nbDofs_uSol + nbDofsPerNode*elementNodes[i]+comp, nbDofs_uSol + nbDofs_uFlu + elementNodes[j], -1.0*Bij, ADD_VALUES);
          }
          if (nbDofs_uSol + nbDofs_uFlu + elementNodes[j] >=m && nbDofs_uSol + nbDofs_uFlu + elementNodes[j] < n){
            double Bij = fe.mixed(i,j,comp);
            MatSetValue(A, nbDofs_uSol + nbDofs_uFlu + elementNodes[j], nbDofsPerNode * elementNodes[i] + comp, par.timeStep() * Bij, ADD_VALUES);
            MatSetValue(A, nbDofs_uSol + nbDofs_uFlu + elementNodes[j], nbDofs_uSol + nbDofsPerNode * elementNodes[i] + comp, Bij, ADD_VALUES);
          }
        }
        if (nbDofs_uSol + nbDofs_uFlu + elementNodes[i] >= m && nbDofs_uSol + nbDofs_uFlu + elementNodes[i] < n){
          double Kij = fe.stiffness(i,j);
          double Mij = fe.mass(i,j);
          MatSetValue(M, nbDofs_uSol + nbDofs_uFlu + elementNodes[i], elementNodes[j], Mij, ADD_VALUES);
          MatSetValue(A, nbDofs_uSol + nbDofs_uFlu + elementNodes[i], nbDofs_uSol + nbDofs_uFlu + elementNodes[j], fe.size()*fe.size() * Kij, ADD_VALUES);
        }
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
    if (world_rank == 0) cout << "Solving BIOT equations for time: " << t << endl;
    if (world_rank == 0) cout << "Iteration: " << t/par.timeStep() << endl;

    /* compute div(u0) */
    Vec uSol_vec = zeros(nbDofs_uSol);
    int low, high;
    VecGetOwnershipRange(uSol_vec, &low, &high);
    Vec u0_seq = getSequential(u0);
    for (int i = 0; i < nbDofs_uSol; i++){
      double uii;
      if (i >= low && i < high){
        code = VecGetValues(u0_seq, 1, &i, &uii); CHKERRQ(code);
        code = VecSetValue(uSol_vec, i, uii, INSERT_VALUES); CHKERRQ(code);
      }
    }
    code = VecAssemblyBegin(uSol_vec); CHKERRQ(code);
    code = VecAssemblyEnd(uSol_vec); CHKERRQ(code);

    /* 2D 3D */
    Vec div_u0 = zeros(geo.nbVertices);
    div_u0 = cal.divergence(uSol_vec);
    double * divU0arr = new double[nbDofsPress];
    Vec b0 = zeros(nbDofs);
    Vec div_u0_seq = getSequential(div_u0);
    code = VecGetValues(div_u0_seq, nbDofsPress, &range(nbDofsPress)[0], divU0arr); CHKERRQ(code);
    for (int i = 0; i < nbDofsPress; i++){
      if (i >= m && i < n){
        code = VecSetValue(b0, i, divU0arr[i], INSERT_VALUES); CHKERRQ(code);
      }
    }
    code = VecAssemblyBegin(b0); CHKERRQ(code);
    code = VecAssemblyEnd(b0); CHKERRQ(code);
    double norm_div_u0;
    VecNorm(div_u0, NORM_2, &norm_div_u0);
    if (world_rank == 0) cout << "BIOT: norm(div u0) = " << norm_div_u0 << endl;

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

    io.writeState(u, t + par.timeStep());

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
