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

#include<heat.hpp>

void Heat::initialize(Parameters parameters, const Geometry & geometry, const Boundary & boundary){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "Heat: Initializing" << endl;
  par = parameters;
  geo = geometry;
  bd = boundary;
  fe.initialize(par, geo.dimension());
  m_verbose = par.verbose();
  nbDofs = 0;
  assert(par.nbDofsPerNode().size() == 1);

  nbDofsPerNode = par.nbDofsPerNode()[0]; 
  nbVertices = geo.nbVertices;
  nbDofs = nbVertices * nbDofsPerNode;

  u0 = zeros(nbDofs);

//  nbNodesPerElement = geo.dimension()+1; 
  mat(A, nbDofs, nbDofs); 
  mat(M, nbDofs, nbDofs); 
//  code = MatGetOwnershipRange(A, &m, &n); CHKERR(code);

  femStat.initialize(par, geo, bd, A);
  femMass.initialize(par, geo, bd, M);
}

void Heat::finalize(){
  MatDestroy(&M);
  VecDestroy(&u0);
  KSPDestroy(&ksp);
}

Mat Heat::assembleLHS(){

//  /* local values */
//  int nbNodesPerElement = geo.dimension()+1; 
//
//  int m,n;
//  code = MatGetOwnershipRange(A, &m, &n); CHKERR(code);

  if (m_world_rank == 0) cout << "Laplacian: Assembling discretization matrix." << endl;
  for (int partId = 0; partId < geo.elements().size(); partId++){
    for (vector<int> simplex : geo.elements()[partId]){ 
      femStat.setSimplex(simplex);
      femMass.copySimplex(femStat);
      for (int var = 0; var < par.nbVariables(); var++){
        femMass.u_dot_v_scalar(1.0, var, var);

        for (int comp = 0; comp < geo.dimension(); comp++){
          femStat.du_x_dv_y(par.viscosity(), comp, comp, var, var);
        }
      }
    }
  }
//
//  for (int partId = 0; partId < geo.elements().size(); partId++){
//
//    if (m_world_rank == 0) cout << "Heat: Assembling discretization matrix. Part " << partId << endl;
//
//    MasterElement fe; 
//    fe.initialize(par, geo.dimension());
//
//    for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on finiteElement */
//      
//      /* get finite element coordinates */
//      vector<int> simplex = geo.elements()[partId][feId];
//      vector<vector<double>> coordinates(nbNodesPerElement);
//      for (int nodeId = 0; nodeId < nbNodesPerElement; nodeId++){
//        coordinates[nodeId] = geo.coordinates()[simplex[nodeId]];
//      }
//      fe.setCoordinates(coordinates);
//
//      if (feId % (geo.elements()[partId].size() / par.verbose()) == 0){
//        if (m_world_rank == 0) cout << "    Elementary matrix for element: " << feId << "/" << geo.elements()[partId].size() - 1 << endl;
//        if (m_world_rank == 0) cout << "    detJac: " << fe.detJacobian() << endl;
//      }
//      /* Assemble elementary matrices */
//      for (int i = 0; i < nbNodesPerElement; i++){
//        for (int j = 0; j < nbNodesPerElement; j++){
//
//          /* u \cdot v */
//          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
//            if (nbDofsPerNode*simplex[i]+comp >= m && nbDofsPerNode*simplex[i]+comp < n){
//              double Mij = par.density() * fe.mass(i,j);
//              code = MatSetValue(M, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp, Mij, ADD_VALUES); CHKERR(code);
//            }
//          }
//
//          /* \grad(u) : \grad(v)  */
//          for (int comp = 0; comp < nbDofsPerNode; comp++){
//            if (nbDofsPerNode*simplex[i]+comp >= m && nbDofsPerNode*simplex[i]+comp < n){
//              double Kij = fe.stiff(i,j); /* diagonal components */
//              code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp,  par.viscosity()*Kij, ADD_VALUES); CHKERR(code);
//            }
//          }
//        }
//      }
//    }
//  }

  code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERR(code);

  /* LHS static matrix */
  Mat K = mat(nbDofs, nbDofs);
  code = MatDuplicate(A, MAT_COPY_VALUES, &K); CHKERR(code);
  code = MatScale(M, 1.0/(par.timeStep())); CHKERR(code);
  code = MatAXPY(K, 1.0, M, DIFFERENT_NONZERO_PATTERN); CHKERR(code);

  double normLHS = norm(K);
  if (m_world_rank == 0) cout << "Heat: Norm LHS = " << normLHS << endl;

  return K;
}

Vec Heat::assembleRHS(){
  Vec b = zeros(nbDofs);
  code = MatMult(M, u0, b); CHKERR(code);
  double normRHS  = norm(b);
  if (m_world_rank == 0) cout << "Heat: Norm RHS = " << normRHS << endl;
  return b;
}

void Heat::update(Vec u, double time){
  code = VecZeroEntries(u0); CHKERR(code);
  code = VecAXPY(u0, 1.0, u); CHKERR(code);
}

void Heat::setSolver(){
  if (m_world_rank == 0) cout << "KSP: Configuring" << endl;
  configureKSP(ksp, par);
}

void Heat::setLHS(Mat C){
  if (m_world_rank == 0) cout << "KSP: Seting operators" << endl;
  code = KSPSetOperators(ksp, C, C); CHKERR(code); 
  code = KSPSetUp(ksp); CHKERR(code);
}

void Heat::solve(Vec b, Vec u){
  if (m_world_rank == 0) cout << "Heat: Solving linear system" << endl;
  code = KSPSolve(ksp, b, u); CHKERR(code); 
  int its;
  double rnorm;
  code = KSPGetResidualNorm(ksp, &rnorm); CHKERR(code);
  code = KSPGetIterationNumber(ksp, &its); CHKERR(code);
  if (par.solver() == "gmres"){
    if (m_world_rank == 0) cout << "Heat: " << its << " iterations. Residual norm = " << rnorm << endl;
  }
}
