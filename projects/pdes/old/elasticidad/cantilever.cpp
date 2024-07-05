/*=============================================================================
  This file is part of the code MAD
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2021,
    
     Felipe Galarce

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
vector<double> stress(vector<double> x, double t){
  vector<double> bc(3, 0.0);
  bc[0] = -1.0;
//  bc[0] = -x[2];
  return bc;
}

vector<double> empotrado(vector<double> x, double t){
  vector<double> bc(3, 0.0);
  return bc;
}

int main(int argc, char *argv[]){

  /* Deploy petsc and slepc */
  PetscInitialize(&argc, &argv, NULL, NULL); 
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

  int nbDofsPerNode = par.nbDofsPerNode()[0]; 
  int nbNodesPerElement = io.dimension()+1; 
  int nbDofs = io.nbVertices() * nbDofsPerNode;

  Mat A = mat(nbDofs, nbDofs); /* Static matrix */
  Vec u = zeros(nbDofs);

  /* Boundary conditions */
  Boundary boundary;
  boundary.initialize(par, geo);
  boundary.Dirichlet(1, empotrado);
//  boundary.Neumann(2, stress);
  boundary.Dirichlet(6, stress);

  int m,n;
  code = MatGetOwnershipRange(A, &m, &n); CHKERRQ(code);

  /* Compute Lame coefficients */
  double lambda = par.youngModulus()/(2.0*(1.0+par.poissonRatio()));
  double mu = par.youngModulus() * par.poissonRatio() / ((1.0 + par.poissonRatio())*(1.0 - 2.0*par.poissonRatio()));

  if (world_rank == 0) cout << "Elasticity: Lame parameters. mu = " << mu << ". lambda = " << lambda << ". " << endl;

  /* 2D 3D */
  MasterTetrahedron finiteElement; 
  finiteElement.initialize(par);

  if (world_rank == 0) cout << "Elasticity: Assembling discretization matrix" << endl;

  for (int finiteElementId = 0; finiteElementId < geo.tetrahedron()[0].size(); finiteElementId++){ /* loop on finiteElement */
    
    /* get finite element coordinates */
    vector<int> elementNodes = geo.tetrahedron()[0][finiteElementId];
    vector<double> mean_point(geo.dimension());
    vector<vector<double>> coordinates(nbNodesPerElement);
    for (int nodeId = 0; nodeId < nbNodesPerElement; nodeId++){
      coordinates[nodeId] = geo.coordinates()[elementNodes[nodeId]];
    }
    finiteElement.setCoordinates(coordinates);
    finiteElement.computeSize();

    if (finiteElementId % (geo.tetrahedron()[0].size() / 5) == 0){
      if (world_rank == 0) cout << "  Elementary matrix for tetra: " << finiteElementId << "/" << geo.tetrahedron()[0].size() - 1 << endl;
      if (world_rank == 0) cout << "  Size: " << finiteElement.size() << endl;
      if (world_rank == 0) cout << "  detJac: " << finiteElement.detJacobian() << endl;
    }
    /* Assemble elementary matrices */
    for (int i = 0; i < nbNodesPerElement; i++){
      for (int j = 0; j < nbNodesPerElement; j++){
        for (int comp = 0; comp < nbDofsPerNode; comp++){
          if (nbDofsPerNode*elementNodes[i]+comp >= m && nbDofsPerNode*elementNodes[i]+comp < n){
            double Kij = finiteElement.stiffness_symmetric(i,j);
            double div_div_ij = finiteElement.div_div(i,j,comp);
            code = MatSetValue(A, nbDofsPerNode*elementNodes[i]+comp, nbDofsPerNode*elementNodes[j]+comp, 
                  mu*Kij + lambda*div_div_ij, 
                  ADD_VALUES); CHKERRQ(code);
          }
        }
      }
    }
  }
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

  /* block and assembly */
  Vec b = zeros(nbDofs);
  boundary.block(b);
  boundary.block(A);

  /* Krylov solver */
  KSP ksp;
  configureKSP(ksp, A, "preonly", "lu");

  if (world_rank == 0) cout << "\n - - - - - - - - - - - - - - - - - - - " << endl;
  if (world_rank == 0) cout << "Solving Elasticity equations. " << endl;
  code = KSPSolve(ksp, b, u); CHKERRQ(code); 

  io.writeState(u);
//  io.writeState(b);
  double norm_sol = norm(u);
  if (world_rank == 0) cout << "Elasticity: Norm solution = " << norm_sol << endl;

  KSPDestroy(&ksp);
  MatDestroy(&A);
  VecDestroy(&u);
  par.finalize();
  PetscFinalize();
}
