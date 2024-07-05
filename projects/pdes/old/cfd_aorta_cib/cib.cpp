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

/* function pointers are passed as argument to BC class */
double inlet_shape(vector<double> x, double t, vector<double> center, double r){
  return 1.0*(1.0 - dot(center - x, center - x) / (r*r));
}

vector<double> inlet(vector<double> x, double t){
  vector<double> bc(3,0.0);
  bc[2] = 1.0;
  return bc;
}

vector<double> noslip(vector<double> x, double t){
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

  CFD cfd;
  cfd.initialize(par, geo, io);

  int nbDofs = (par.nbDofsPerNode()[0] + par.nbDofsPerNode()[1])*io.nbVertices();
  int nbDofsPress = io.nbVertices();
  int nbDofsVel = par.nbDofsPerNode()[0]*io.nbVertices();
  int nbVertices = io.nbVertices();
  int nbDofsPerNode = par.nbDofsPerNode()[0];

  /* Save initial condition */
  Vec phi0 = zeros(nbDofs);
  Vec phi = zeros(nbDofs);
  io.writeState(phi0, 0.0);

  int m_vel, n_vel;
  VecGetOwnershipRange(boundary.normals(), &m_vel, &n_vel);

  Mat M = mat(nbDofs, nbDofs); /* Mass matrix */
  Mat B = mat(nbDofs, nbDofs); /* Static matrix */

  int m,n;
  code = MatGetOwnershipRange(M, &m, &n); CHKERRQ(code);
  if (world_rank == 0) cout << "NS: Assembling discretization matrix." << endl;
  for (int tetraId = 0; tetraId < geo.tetrahedron()[0].size(); tetraId++){ /* loop on tetra */

    vector<int> tetra = geo.tetrahedron()[0][tetraId];
    /* get finite element coordinates */
    vector<vector<double>> coordinates(4);
    coordinates[0] = geo.coordinates()[tetra[0]];
    coordinates[1] = geo.coordinates()[tetra[1]];
    coordinates[2] = geo.coordinates()[tetra[2]];
    coordinates[3] = geo.coordinates()[tetra[3]];

    tet.setCoordinates(coordinates);
    tet.computeSize();

    if (tetraId % (geo.tetrahedron()[0].size() / 5) == 0){
      if (world_rank == 0) cout << "  Elementary matrix for tetra: " << tetraId << "/" << geo.tetrahedron()[0].size() - 1 << endl;
    }
    /* Assemble elementary matrices */
    for (int j = 0; j < 4; j++){
      for (int i = 0; i < 4; i++){
        for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
          if (par.nbDofsPerNode()[0]*tetra[i]+comp >= m && par.nbDofsPerNode()[0]*tetra[i]+comp < n){
            double Kij = tet.stiffness_symmetric(i,j);
            double Mij = par.density()/par.timeStep()*tet.mass(i,j);
            /* stiffness */
            MatSetValue(M, par.nbDofsPerNode()[0]*tetra[i]+comp, par.nbDofsPerNode()[0]*tetra[j]+comp, Mij, ADD_VALUES); 
            MatSetValue(B, par.nbDofsPerNode()[0]*tetra[i]+comp, par.nbDofsPerNode()[0]*tetra[j]+comp, par.viscosity() * Kij, ADD_VALUES); 
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
          double Kij = tet.stiffness_symmetric(i,j);
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

  /* LHS static matrix */
  Mat A = mat(nbDofs, nbDofs);
  code = MatDuplicate(M, MAT_COPY_VALUES, &A); CHKERRQ(code);
  code = MatAXPY(A, 1.0, B, DIFFERENT_NONZERO_PATTERN); CHKERRQ(code);

  /* Krylov solver */
  KSP ksp;
  code = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERR(code); 
  PetscBool reuse_prec = PETSC_TRUE;

  double p_out;
  double p_d0 = 0.0;

  for (double t : range(0, par.nbIterations()*par.timeStep(), par.timeStep())){

    if (world_rank == 0) cout << "\n- - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
    if (world_rank == 0) cout << "NS: Solving Navier-Stokes equation for time: " << t << endl;

    /* block and assembly rhs */
    Vec rhs = vec(nbDofs);
    MatMult(M, phi0, rhs);

    /* Apply tension boundary conditions with Windkessel model */
    if (world_rank == 0) cout << "NS: Coupling with 0D solver." << endl;
    Vec tension_out = zeros(nbDofs);

    for (int i = 0; i < par.outlets().size(); i++){
      double Q_out = cfd.flow(phi, par.outlets()[i]);
      double p_d = 1.0/(par.capacitances()[i]/par.timeStep() + 1.0/par.distalResistances()[i]) * (Q_out + p_d0*par.capacitances()[i]/par.timeStep());
      p_out = p_d + par.resistances()[i] * Q_out;
      p_d0 = p_d;

      if (world_rank == 0) cout << "  Boundary label: " << par.outlets()[i] << endl;
      if (world_rank == 0) cout << "    Flow = " << Q_out << "." << endl;
      if (world_rank == 0) cout << "    Distal pressure = " << p_d << "." << endl;
      if (world_rank == 0) cout << "    Proximal pressure = " << par.resistances()[i] * Q_out << "." << endl;

      for (int ind = 0; ind < nbDofsVel; ind++){
        double n_ij;
        if (ind >= m_vel && ind < n_vel){
          code = VecGetValues(boundary.normals(), 1, &ind, &n_ij); CHKERRQ(code);
          code = VecSetValue(tension_out, ind, 1.0*p_out*n_ij, ADD_VALUES); CHKERRQ(code);
        }
      }
    }
    code = VecAssemblyBegin(tension_out); CHKERR(code);
    code = VecAssemblyEnd(tension_out); CHKERR(code);
    double normTension;
    VecNorm(tension_out, NORM_2, &normTension);
    if (world_rank == 0) cout << "    Norm of boundary tension = " << normTension << "." << endl;
    Vec rhs_windkessel = vec(nbDofs);
    code = MatMult(boundary.mass(), tension_out, rhs_windkessel); CHKERR(code);
    code = VecAXPY(rhs, 1.0, rhs_windkessel); CHKERR(code);

    /* Assemble time dependent matrix */
    Mat C = mat(nbDofs, nbDofs);

    Vec phi0_seq = getSequential(phi);

    if (world_rank == 0) cout << "NS: Assembling convection matrix." << endl;
    for (int tetraId = 0; tetraId < geo.tetrahedron()[0].size(); tetraId++){ /* loop on tetra */

      vector<int> tetra = geo.tetrahedron()[0][tetraId];
      /* get finite element coordinates */
      vector<vector<double>> coordinates(4);
      coordinates[0] = geo.coordinates()[tetra[0]];
      coordinates[1] = geo.coordinates()[tetra[1]];
      coordinates[2] = geo.coordinates()[tetra[2]];
      coordinates[3] = geo.coordinates()[tetra[3]];

      tet.setCoordinates(coordinates);
      tet.computeSize();

      if (tetraId % (geo.tetrahedron()[0].size() / 5) == 0){
        if (world_rank == 0) cout << "  Elementary matrix for tetra: " << tetraId << "/" << geo.tetrahedron()[0].size() - 1 << endl;
      }

      /* loc2glob dof mapping */
      vector<int> loc2glob(12);
      for (int i = 0; i < 4; i++){
        for (int comp = 0; comp < 3; comp++){
          loc2glob[3*i+comp] = 3*tetra[i]+comp;
        }
      }

      /* Gather advection field on element */
      vector<double> u_el(12);
      code = VecGetValues(phi0_seq, 12, &loc2glob[0], &u_el[0]); CHKERR(code);

      /* Set values from local to global */
      for (int j = 0; j < 4; j++){
        vector<double> u_node_j(par.nbDofsPerNode()[0]);
        for (int i = 0; i < 4; i++){
          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
            u_node_j[comp] = u_el[par.nbDofsPerNode()[0]*j+comp];
          }
          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
            if (par.nbDofsPerNode()[0]*tetra[i]+comp >= m && par.nbDofsPerNode()[0]*tetra[i]+comp < n){
              double Aij = par.density()*tet.advection(i, j, u_node_j);
              MatSetValue(C, par.nbDofsPerNode()[0]*tetra[i]+comp, par.nbDofsPerNode()[0]*tetra[j]+comp, Aij, ADD_VALUES);
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
    //boundary.Dirichlet(par.inlet(), inlet, 0);
    boundary.DirichletNormal(par.inlet(), inlet_shape, 0);
    for (int i = 0; i < par.walls().size(); i++){
      boundary.Dirichlet(par.walls()[i], noslip, 0);
    }
    boundary.block(rhs);
    boundary.block(C);

    code = KSPSetOperators(ksp, C, C); CHKERR(code); 
    code = KSPSetType(ksp, KSPPREONLY); CHKERR(code); 
    PC pc;
    code = KSPGetPC(ksp, &pc); CHKERR(code); 
    code = PCSetType(pc, PCLU); CHKERR(code); 
    code = KSPSetReusePreconditioner(ksp, reuse_prec); CHKERR(code);
    code = KSPSetFromOptions(ksp); CHKERR(code); 

    if (world_rank == 0) cout << "KSP: Solving linear system" << endl;
    code = KSPSolve(ksp, rhs, phi); CHKERRQ(code); 
    int its;
    double rnorm;
    code = KSPGetResidualNorm(ksp,&rnorm); CHKERRQ(code);
    code = KSPGetIterationNumber(ksp,&its); CHKERRQ(code);
    if (world_rank == 0) cout << "KSP: " << its << " iterations. Residual norm = " << rnorm << endl;
  
    for (int idBD : par.outlets()){
      boundary.backflow(phi);
    }

    io.writeState(phi, t + par.timeStep());

    VecZeroEntries(phi0);
    code = VecAXPY(phi0, 1.0, phi); CHKERRQ(code);

    double norm_sol = norm(phi);
    if (world_rank == 0) cout << "KSP: Norm solution = " << norm_sol << endl;

    VecDestroy(&phi0_seq);
    VecDestroy(&rhs);
    VecDestroy(&rhs_windkessel);
    VecDestroy(&tension_out);
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
