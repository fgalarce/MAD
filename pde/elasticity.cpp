/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2021,
    
     Felipe Galarce at INRIA

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

#include<elasticity.hpp>

void Elasticity::initialize(Parameters parameters, const Geometry & geometry){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "Elasticity: Initializing" << endl;
  par = parameters;
  geo = geometry;
  fe.initialize(par, geo.dimension());
  m_verbose = par.verbose();
  nbDofs = 0;
  assert(par.nbDofsPerNode().size() == 1);

  nbDofsPerNode = par.nbDofsPerNode()[0]; 
  nbVertices = geo.nbVertices;
  nbDofs = nbVertices * nbDofsPerNode;

}

void Elasticity::finalize(){
  MatDestroy(&M);
  KSPDestroy(&ksp);
}

Mat Elasticity::assembleLHS(){

  mat(A, nbDofs, nbDofs); 
  mat(M, nbDofs, nbDofs); 

  /* local values */
  int nbNodesPerElement = geo.dimension()+1; 

  int m,n;
  code = MatGetOwnershipRange(A, &m, &n); CHKERR(code);

  for (int partId = 0; partId < geo.elements().size(); partId++){

    if (m_world_rank == 0) cout << "Elasticity: Assembling discretization matrix. Part " << partId << endl;
      /* Compute Lame coefficients */
    double lambda = par.youngModulus()[partId]/(2.0*(1.0+par.poissonRatio()[partId]));
    double mu = par.youngModulus()[partId] * par.poissonRatio()[partId] / ((1.0 + par.poissonRatio()[partId])*(1.0 - 2.0*par.poissonRatio()[partId]));
    if (m_world_rank == 0) cout << "Elasticity: Lame parameters. mu = " << mu << ". lambda = " << lambda << ". " << endl;

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
        if (m_world_rank == 0) cout << "    Elementary matrix for element: " << feId << "/" << geo.elements()[partId].size() - 1 << endl;
        if (m_world_rank == 0) cout << "    detJac: " << fe.detJacobian() << endl;
      }
      /* Assemble elementary matrices */
      for (int i = 0; i < nbNodesPerElement; i++){
        for (int j = 0; j < nbNodesPerElement; j++){

          /* u \cdot v */
          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
            if (nbDofsPerNode*simplex[i]+comp >= m && nbDofsPerNode*simplex[i]+comp < n){
              double Mij = par.density() * fe.mass(i,j);
              code = MatSetValue(M, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp, Mij, ADD_VALUES); CHKERR(code);
            }
          }

          /* \epsilon(u) : \epsilon(v)  */
          for (int comp = 0; comp < nbDofsPerNode; comp++){
            if (nbDofsPerNode*simplex[i]+comp >= m && nbDofsPerNode*simplex[i]+comp < n){
              double Kij = 1.0*fe.stiffness_symmetric(i,j); /* diagonal components */
              code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp, mu*Kij, ADD_VALUES); CHKERR(code);
            }
          }
          for (int comp_v = 0; comp_v < nbDofsPerNode; comp_v++){
            if (nbDofsPerNode*simplex[i]+comp_v>= m && nbDofsPerNode*simplex[i]+comp_v < n){
              for (int comp_u = 0; comp_u < nbDofsPerNode; comp_u++){
                double Kdev = 1.0/2.0 * fe.dphi_dx(i, comp_u) * fe.dphi_dx(j, comp_v) * fe.detJacobian()*fe.weightsSum();
                code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp_v, nbDofsPerNode*simplex[j]+comp_u, mu*Kdev, ADD_VALUES); CHKERR(code);
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

  code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERR(code);

  /* LHS static matrix */
  double normLHS = norm(A);
  if (m_world_rank == 0) cout << "Elasticity: Norm LHS = " << normLHS << endl;

  return A;
}

Vec Elasticity::assembleRHS(){
  Vec b = zeros(nbDofs);
  /* TODO: implement gravity forces */
  double normRHS  = norm(b);
  if (m_world_rank == 0) cout << "Elasticity: Norm RHS = " << normRHS << endl;
  return b;
}

void Elasticity::setSolver(){
  if (m_world_rank == 0) cout << "KSP: Configuring" << endl;
  configureKSP(ksp, par);
}

void Elasticity::setLHS(Mat C){
  if (m_world_rank == 0) cout << "KSP: Seting operators" << endl;
  code = KSPSetOperators(ksp, C, C); CHKERR(code); 
  code = KSPSetUp(ksp); CHKERR(code);
}

void Elasticity::solve(Vec b, Vec u){
  if (m_world_rank == 0) cout << "Elasticity: Solving linear system" << endl;
  code = KSPSolve(ksp, b, u); CHKERR(code); 
  int its;
  double rnorm;
  code = KSPGetResidualNorm(ksp, &rnorm); CHKERR(code);
  code = KSPGetIterationNumber(ksp, &its); CHKERR(code);
  if (par.solver() == "gmres"){
    if (m_world_rank == 0) cout << "Elasticity: " << its << " iterations. Residual norm = " << rnorm << endl;
  }
}
