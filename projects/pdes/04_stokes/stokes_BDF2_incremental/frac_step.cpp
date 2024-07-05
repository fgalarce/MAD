/*=============================================================================
  This file is part of the code M.D.M.A. (or MDMA)
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2020,
    
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

#include <ultra-4d-flow.hpp>

vector<double> noslip(vector<double> x, double t){
  vector<double> bc(3,0.0);
  return bc;
}

vector<double> inlet(vector<double> x, double t){
  vector<double> bc(3, 0.0);
//  bc[1] = (1 - (x[0]*x[0] + x[2]*x[2])/(0.2*0.2))*sin(2*PI*t);
  bc[1] = 1.0;
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

  /* Initialize MDMA objects */
  IO io;
  io.initialize(par);
    
  Geometry geo;
  geo.initialize(par, io);

  CFD cfd;
  cfd.initialize(par, geo, io);

  MasterTetrahedron tet; 
  tet.initialize(par);

  BoundaryConditions bc;
  bc.initialize(par, geo);

  int nbDofsPress = io.nbVertices();
  int nbDofsVel = 3*io.nbVertices();

  Mat K_p = mat(nbDofsPress, nbDofsPress); /* Scalar stiffness */
  Mat M_p = mat(nbDofsPress, nbDofsPress); /* Scalar mass */
  Mat M = mat(nbDofsVel, nbDofsVel); /* Mass */
  Mat K = mat(nbDofsVel, nbDofsVel); /* Stiffness */

  int m,n;
  code = MatGetOwnershipRange(M, &m, &n); CHKERR(code);

  Tic tic;
  int tetraId = 0; 
  for (vector<int> tetra : geo.tetrahedron()[0]){ /* Loop on tetra */
    /* get finite element coordinates */
    vector<vector<double>> coordinates(4);
    coordinates[0] = geo.coordinates()[tetra[0]];
    coordinates[1] = geo.coordinates()[tetra[1]];
    coordinates[2] = geo.coordinates()[tetra[2]];
    coordinates[3] = geo.coordinates()[tetra[3]];

    tet.setCoordinates(coordinates);

    if (tetraId % (geo.tetrahedron()[0].size() / 10) == 0){
      if (world_rank == 0) cout << "  Computing elementary arrays for tetrahedron " << tetraId << " of " << geo.tetrahedron()[0].size() - 1 << endl;
      if (world_rank == 0) cout << "  Labels: " << tetra << endl; 
      if (world_rank == 0) cout << "  det(J) = " << tet.detJacobian() << endl; 
    }
    /* Assemble elementary matrices */
    for (int i = 0; i < 4; i++){
      for (int j = 0; j < 4; j++){

        if (tetra[i] >= m && tetra[i] < n){
          code = MatSetValue(M_p, tetra[i], tetra[j], tet.mass(i,j), ADD_VALUES); CHKERRQ(code);
          code = MatSetValue(K_p, tetra[i], tetra[j], tet.stiffness(i,j), ADD_VALUES); CHKERRQ(code);
        }
        for (int comp = 0; comp < 3; comp++){
          if ((3*tetra[i]+comp) >= m && (3*tetra[i]+comp) < n){
            code = MatSetValue(M, 3*tetra[i]+comp, 3*tetra[j]+comp, tet.mass(i,j), ADD_VALUES); CHKERRQ(code);
            code = MatSetValue(K, 3*tetra[i]+comp, 3*tetra[j]+comp, tet.stiffness(i,j), ADD_VALUES); CHKERRQ(code);
          }
        }
      }
    }
    tetraId++;
  }
  if (world_rank == 0) cout << " Time elapsed: " << tic.toc() << " sec. "<< endl;
  code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyBegin(K_p, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(K_p, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyBegin(M_p, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(M_p, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

  /* Velocity */
  Vec u3 = vec(nbDofsVel);
  Vec u0 = vec(nbDofsVel); 
  Vec u1 = vec(nbDofsVel); /* guess step */
  Vec u2 = vec(nbDofsVel); /* update step */

  /* Pressure */
  Vec p0 = zeros(nbDofsPress);
  Vec p1 = vec(nbDofsPress);

  /* save initial condition */
//  io.writeState(u0, "Gr_p", 0.0);
//  io.writeState(u1, "velocity", 0.0);
//  io.writeState(u0, "velocity0", 0.0);
//  io.writeState(p1, "pressure", 0.0);
//  io.writeState(p1, "divu", 0.0);

  for (double t : rango(par.timeStep(), par.nbIterations()*par.timeStep(), par.timeStep())){
    if (world_rank == 0) cout << "\n- - - - - - - - - - - - - - - - - - - - - - - - " << endl;
    if (world_rank == 0) cout << "Solving stokes equation for time: " << t << endl;

    if (world_rank == 0) cout << "  Viscous step A u = a. " << endl;
    /* Linear system: A u = a */
    /*                (3/2 * M + dt * nu * K) u^* = 2 * M * u0 - 1/2 * M * u3 - dt * M \grad p0 */
    Mat A = mat(nbDofsVel, nbDofsVel);
    Vec a = vec(nbDofsVel);
    Vec a0 = vec(nbDofsVel);
    Vec a1 = vec(nbDofsVel);

    /* Prepare A */
    code = MatDuplicate(M, MAT_COPY_VALUES, &A); CHKERRQ(code);
    code = MatScale(A, 3.0/2.0); CHKERRQ(code);
    code = MatAXPY(A, par.timeStep() * par.viscosity(), K, SAME_NONZERO_PATTERN); CHKERRQ(code);
  
    /* Prepare a */
    code = MatMult(M, u0, a); CHKERRQ(code);
    code = VecScale(a, 2.0); CHKERRQ(code);
    code = MatMult(M, cfd.gradient(p0), a0); CHKERRQ(code);
    code = VecAXPY(a, -1.0*par.timeStep(), a0); CHKERRQ(code);
    code = MatMult(M, u3, a1); CHKERRQ(code);
    code = VecAXPY(a, -1.0/2.0, a1); CHKERRQ(code);
    
    
    /* Apply boundary conditions */
    if (world_rank == 0) cout << "    Norm(a) " << norm(a) << endl;
    if (world_rank == 0) cout << "    Norm(A) " << norm(A) << endl;
    if (world_rank == 0) cout << "    Applying boundary conditions." << endl;
    bc.time(t);
    bc.Dirichlet(par.inlet(), inlet);
//    bc.Dirichlet(par.inlet(), "logit_inlet");
    bc.Dirichlet(par.wall(), noslip);

    /* block and assembly */
    bc.block(A);
    bc.block(a); 
    if (world_rank == 0) cout << "    l2-norm of a: " << norm(a) << endl;
    if (world_rank == 0) cout << "    Fro-norm of A: " << norm(A) << endl;

    /* Solve */
    KSP solver = configureKSP(A, "preonly", "lu");
    code = KSPSolve(solver, a, u1); CHKERRQ(code); 

//    io.writeState(u1, "velocity0", t - par.timeStep());

    /* clean */
    KSPDestroy(&solver);
    VecDestroy(&a);
    MatDestroy(&A);

    if (world_rank == 0) cout << "Projection into Hdiv. Pressure. " << endl;
    /* Linear System B p = b */
    /*              dt * K_p * p = dt * K_p * p0 - M_p * div(u1) */
    Mat B = mat(nbDofsPress, nbDofsPress);
    Vec b = zeros(nbDofsPress);
    Vec b0 = zeros(nbDofsPress);

    /* Prepare b */
    Vec divu = cfd.divergence(u1);
//    io.writeState(divu, "divu", t - par.timeStep());

    double norm_divu = sqrt(ip(divu, M_p, divu));
    if (world_rank == 0) cout << "  L2 norm of hat u divergence: " << norm_divu << endl;
    code = MatMult(M_p, divu, b0); CHKERRQ(code);
    code = MatMult(K_p, p0, b); CHKERRQ(code);
    code = VecScale(b, par.timeStep());
    code = VecAXPY(b, -3.0/2.0, b0);

    /* Prepare B */
    code = MatDuplicate(K_p, MAT_COPY_VALUES, &B); CHKERRQ(code);
    code = MatScale(B, par.timeStep());
    /* Outlet bc for pressure */
    int bdLabel = par.outlets()[0];
    int nbDofsBC = geo.boundaryNodes(bdLabel).size();
    vector<double> BC(nbDofsBC, 0.0);
    vector<int> iBC = geo.boundaryNodes(bdLabel);
    code = VecSetValues(b, nbDofsBC, &iBC[0], &BC[0], INSERT_VALUES); CHKERRQ(code);
    code = MatZeroRows(B, nbDofsBC, &iBC[0], 1.0, NULL, NULL); CHKERRQ(code);

    /* Assembly */
    code = VecAssemblyBegin(b); CHKERR(code); CHKERRQ(code);
    code = VecAssemblyEnd(b); CHKERR(code); CHKERRQ(code);
    code = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
    code = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

    /* Solve */
    solver = configureKSP(B, "preonly", "lu");
    code = KSPSolve(solver, b, p1); CHKERRQ(code); 

    /* clean */
    KSPDestroy(&solver);
    VecDestroy(&b);
    MatDestroy(&B);

    io.writeState(p1, "pressure", t - par.timeStep());

    if (world_rank == 0) cout << "Projection into Hdiv. Velocity. " << endl;
    Vec grad_p = cfd.gradient(p1);
//    io.writeState(grad_p, "Gr_p", t - par.timeStep());
    VecCopy(u1, u2);
    VecAXPY(u2, -1.0*par.timeStep(), cfd.gradient(p1));
    VecAXPY(u2,  1.0*par.timeStep(), cfd.gradient(p0));
    bc.block(u2);

    io.writeState(u2, "velocity", t - par.timeStep());
    VecCopy(u0,u3);
    VecCopy(u2,u0);
    VecCopy(p1,p0);
  }

  VecDestroy(&u0);
  VecDestroy(&u1);
  VecDestroy(&u2);
  VecDestroy(&p1);
  MatDestroy(&M);
  MatDestroy(&K);
  MatDestroy(&K_p);
}
