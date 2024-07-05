/*=============================================================================
  This file is part of the code MAD
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2017/2021,
    
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

double MRE(vector<double> x, double t, Parameters par){
  double r = 7.5;
  double bc = par.amplitude() * ( sin(2*PI * t / par.period()) ) * ( 1.0 - ((r - x[1])*(r - x[1])) / (r*r) );
//  double bc = par.amplitude() * ( 1.0 - ((r - x[1])*(r - x[1])) / (r*r) );
//  double bc = par.amplitude() * ( sin(2*PI * t / par.period()) );
  if (bc > 0) {
    bc = 0.0;
  }
  return bc;
}

vector<double> neck(vector<double> x, double t, Parameters par){
  vector<double> bc(2, 0.0);
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
  Mat M = mat(nbDofs, nbDofs); /* Mass matrix */

  Vec u00 = zeros(nbDofs);
  Vec u0 = zeros(nbDofs);
  Vec u = zeros(nbDofs);

  /* Boundary conditions */
  Boundary boundary;
  boundary.initialize(par, geo);

  io.writeState(u0, 0.0);

  int m,n;
  code = MatGetOwnershipRange(A, &m, &n); CHKERRQ(code);


  for (int partId = 0; partId < geo.elements().size(); partId++){

    if (world_rank == 0) cout << "Elasticity: Assembling discretization matrix. Part " << partId << endl;
      /* Compute Lame coefficients */
    double lambda = par.youngModulus()[partId]/(2.0*(1.0+par.poissonRatio()[partId]));
    double mu = par.youngModulus()[partId] * par.poissonRatio()[partId] / ((1.0 + par.poissonRatio()[partId])*(1.0 - 2.0*par.poissonRatio()[partId]));
    if (world_rank == 0) cout << "Elasticity: Lame parameters. mu = " << mu << ". lambda = " << lambda << ". " << endl;

    MasterElement fe; 
    fe.initialize(par, geo.dimension());

    for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on finiteElement */
      
      /* get finite element coordinates */
      vector<int> simplex = geo.elements()[partId][feId];
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int nodeId = 0; nodeId < nbNodesPerElement; nodeId++){
        coordinates[nodeId] = geo.coordinates()[simplex[nodeId]];
      }
      fe.setCoordinates(coordinates);

      if (feId % (geo.elements()[partId].size() / par.verbose()) == 0){
        if (world_rank == 0) cout << "    Elementary matrix for element: " << feId << "/" << geo.elements()[partId].size() - 1 << endl;
        if (world_rank == 0) cout << "    detJac: " << fe.detJacobian() << endl;
      }
      /* Assemble elementary matrices */
      for (int i = 0; i < nbNodesPerElement; i++){
        for (int j = 0; j < nbNodesPerElement; j++){

          /* u \cdot v */
          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
            if (nbDofsPerNode*simplex[i]+comp >= m && nbDofsPerNode*simplex[i]+comp < n){
              double Mij = par.density() * fe.mass(i,j) / (par.timeStep()*par.timeStep());
              code = MatSetValue(M, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp, Mij, ADD_VALUES); CHKERR(code);
            }
          }

          /* \epsilon(u) : \epsilon(v)  */
          for (int comp = 0; comp < nbDofsPerNode; comp++){
            if (nbDofsPerNode*simplex[i]+comp >= m && nbDofsPerNode*simplex[i]+comp < n){
              double Kij = 1.0/2.0*fe.stiff(i,j); /* diagonal components */
              code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp,  2.0*mu*Kij, ADD_VALUES); CHKERR(code);
            }
          }
          for (int comp_v = 0; comp_v < nbDofsPerNode; comp_v++){
            if (nbDofsPerNode*simplex[i]+comp_v>= m && nbDofsPerNode*simplex[i]+comp_v < n){
              for (int comp_u = 0; comp_u < nbDofsPerNode; comp_u++){
                double Kdev = 1.0/2.0 * fe.dphi_dx(i, comp_u) * fe.dphi_dx(j, comp_v) * fe.detJacobian()*fe.weightsSum();
                code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp_v, nbDofsPerNode*simplex[j]+comp_u,  2.0*mu*Kdev, ADD_VALUES); CHKERR(code);
              }
            }
          }

          /* div(u) div(v) */
          for (int comp_v = 0; comp_v < nbDofsPerNode; comp_v++){
            /* deviatoric components */
            if (nbDofsPerNode*simplex[i]+comp_v>= m && nbDofsPerNode*simplex[i]+comp_v < n){
              for (int comp_u = 0; comp_u < nbDofsPerNode; comp_u++){
                double div_div_ij = fe.dphi_dx(i, comp_v) * fe.dphi_dx(j, comp_u) * fe.detJacobian()*fe.weightsSum();
                code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp_v, nbDofsPerNode*simplex[j]+comp_u,  lambda*div_div_ij, ADD_VALUES); CHKERR(code);
              }
            }
          }
        }
      }
    }
  }

  code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

  /* A + rho/dt/dt M */
  code = MatAXPY(A, 1.0, M, DIFFERENT_NONZERO_PATTERN); CHKERRQ(code);

  boundary.NeumannNormal(par.bcNeumann()[0], MRE);
  boundary.Dirichlet(par.fixed_boundary(), neck);
  boundary.block(A);

  /* Krylov solver */
  KSP ksp = configureKSP(A, "preonly", "lu");

  ofstream ksp_output(par.dirResults() + "/ksp.log");

  for (double t : range(0, par.nbIterations()*par.timeStep(), par.timeStep())){

    if (world_rank == 0) cout << "\n - - - - - - - - - - - - - - - - - - - " << endl;
    if (world_rank == 0) cout << "Solving Elasticity equations for time: " << t << endl;
    if (world_rank == 0) cout << "Iteration: " << t/par.timeStep() << endl;

    /* block and assembly b */
    boundary.time(t); 
    code = VecAXPY(u0, -1.0, u00); CHKERRQ(code);
    boundary.NeumannNormal(par.bcNeumann()[0], MRE);
    boundary.Dirichlet(par.fixed_boundary(), neck);

    Vec b = zeros(nbDofs);
    code = MatMult(M, u0, b); CHKERRQ(code);
    double Mu0_norm = norm(b);
    if (world_rank == 0) cout << "Elasticity: Norm b=" << Mu0_norm << endl;
    
    boundary.block(b);

    if (world_rank == 0) cout << "Elasticity: Solving linear system" << endl;
    code = KSPSolve(ksp, b, u); CHKERRQ(code); 
    int its;
    double rnorm;
    code = KSPGetResidualNorm(ksp, &rnorm); CHKERRQ(code);
    code = KSPGetIterationNumber(ksp, &its); CHKERRQ(code);

    if (world_rank == 0) cout << "Elasticity: " << its << " iterations. Residual norm = " << rnorm << endl;
    ksp_output << its << " " << rnorm << endl;
    
    io.writeState(u, t + par.timeStep());
       
    /* Update solutions for centered finite differences */ 
    if (t > 0){
      VecZeroEntries(u00); 
      VecAXPY(u00, 1.0, u0);
    }
    VecZeroEntries(u0); 
    code = VecAXPY(u0, 2.0, u); CHKERRQ(code);

    double norm_sol = norm(u);
    if (world_rank == 0) cout << "Elasticity: Norm solution = " << norm_sol << endl;
  }

  ksp_output.close();

  KSPDestroy(&ksp);
  MatDestroy(&A);
  MatDestroy(&M);
  VecDestroy(&u);
  VecDestroy(&u0);
  VecDestroy(&u00);
  par.finalize();
  PetscFinalize();
}
