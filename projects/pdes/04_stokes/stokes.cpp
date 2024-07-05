/*=============================================================================
  This file is part of the code MAD
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2021,
    
     Felipe Galarce at INRIA/WIAS

  MAD is free software; you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by the Free
  Software Foundation; either version 2.1 of the License, or (at your option)
  any later version.

  MAD is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with MDMA. If not, see http://www.gnu.org/licenses/.
  =============================================================================*/

#include <mad.hpp>

vector<double> inlet_time(vector<double> x, double t, Parameters par){
  vector<double> center(3, 0.0);
  vector<double> bc(3, 0.0);
  double r = 0.2;
  bc[1] = par.inlet_u0()*(-1.0*(1.0 - dot(center - x, center - x) / (r*r)))*sin(2*PI*t);
  return bc;
}

vector<double> inlet_cnst(vector<double> x, double t, Parameters par){
  vector<double> center(3, 0.0);
  vector<double> bc(3, 0.0);
  double r = 0.2;
  bc[1] = par.inlet_u0()*(-1.0*(1.0 - dot(center - x, center - x) / (r*r)));
  return bc;
}

vector<double> noslip(vector<double> x, double t, Parameters par){
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

  Boundary boundary;
  boundary.initialize(par, geo);

  int nbDofs = (par.nbDofsPerNode()[0] + par.nbDofsPerNode()[1])*io.nbVertices();
  int nbDofsPress = io.nbVertices();
  int nbDofsVel = par.nbDofsPerNode()[0]*io.nbVertices();
  int nbVertices = io.nbVertices();
  int nbDofsPerNode = par.nbDofsPerNode()[0];
  int nbVerticesPerElement = io.dimension() + 1;

  /* Save initial condition */
  Vec phi0 = zeros(nbDofs);
  Vec phi = zeros(nbDofs);
  io.writeState(phi0, 0.0);

  Mat M = mat(nbDofs, nbDofs); /* Mass matrix */
  Mat B = mat(nbDofs, nbDofs); /* Static matrix */

  int m,n;
  code = MatGetOwnershipRange(M, &m, &n); CHKERRQ(code);
  if (world_rank == 0) cout << "Stokes: Assembling discretization matrix." << endl;
  for (int tetraId = 0; tetraId < geo.tetrahedron()[0].size(); tetraId++){ /* loop on tetra */

    vector<int> tetra = geo.tetrahedron()[0][tetraId];
    /* get finite element coordinates */
    vector<vector<double>> coordinates(nbVerticesPerElement);
    for (int idVer = 0; idVer < nbVerticesPerElement; idVer++){
      coordinates[idVer] = geo.coordinates()[tetra[idVer]];
    }
    tet.setCoordinates(coordinates);
    tet.computeSize();

    if (tetraId % (geo.tetrahedron()[0].size() / 5) == 0){
      if (world_rank == 0) cout << "  Elementary matrix for tetra: " << tetraId << "/" << geo.tetrahedron()[0].size() - 1 << endl;
    }
    /* Assemble elementary matrices */
    for (int i = 0; i < nbVerticesPerElement; i++){
      for (int j = 0; j < nbVerticesPerElement; j++){
        for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
          if (par.nbDofsPerNode()[0]*tetra[i]+comp >= m && par.nbDofsPerNode()[0]*tetra[i]+comp < n){
            double Mij = par.density()/par.timeStep()*tet.mass(i,j);
            double Kij = 2.0 * par.viscosity() * tet.stiffness_symmetric(i,j);
            /* stiffness */
            MatSetValue(M, par.nbDofsPerNode()[0]*tetra[i]+comp, par.nbDofsPerNode()[0]*tetra[j]+comp, Mij, ADD_VALUES); 
            MatSetValue(B, par.nbDofsPerNode()[0]*tetra[i]+comp, par.nbDofsPerNode()[0]*tetra[j]+comp, Kij, ADD_VALUES); 
            /* div_phi_i phi_j */
            double Bij = tet.mixed(i,j,comp);
            MatSetValue(B, par.nbDofsPerNode()[0]*tetra[i]+comp, nbDofsVel + tetra[j], -1.0*Bij, ADD_VALUES);  /* -B*/
          }
          if (nbDofsVel + tetra[j] >= m && nbDofsVel + tetra[j] < n){
            double Bij = tet.mixed(i,j,comp);
            MatSetValue(B, nbDofsVel + tetra[j], par.nbDofsPerNode()[0]*tetra[i]+comp, 1.0*Bij, ADD_VALUES);  /* B^T */
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

  for (double t : range(0, par.nbIterations()*par.timeStep(), par.timeStep())){

    if (world_rank == 0) cout << "\n- - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
    if (world_rank == 0) cout << "Stokes: Solving Navier-Stokes equation for time: " << t << endl;

    /* block and assembly rhs */
    Vec rhs = vec(nbDofs);
    MatMult(M, phi0, rhs);

    boundary.time(t);
    boundary.Dirichlet(par.inlet(), inlet_time, 0);
    for (int i = 0; i < par.walls().size(); i++){
      boundary.Dirichlet(par.walls()[i], noslip, 0);
    }
    boundary.block(rhs);

    /* LHS static matrix */
    Mat A = mat(nbDofs, nbDofs);
    code = MatDuplicate(M, MAT_COPY_VALUES, &A); CHKERRQ(code);
    code = MatAXPY(A, 1.0, B, DIFFERENT_NONZERO_PATTERN); CHKERRQ(code);
    boundary.block(A);
    KSP ksp = configureKSP(A, "preonly", "lu");

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

    VecDestroy(&rhs);
    KSPDestroy(&ksp);
    MatDestroy(&A);
  }

  VecDestroy(&phi);
  VecDestroy(&phi0);
  MatDestroy(&B);
  MatDestroy(&M);
  SlepcFinalize();
}
