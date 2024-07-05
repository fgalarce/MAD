/*=============================================================================
  This file is part of the code MAD
  Copyright (C) 2017-2023,
    
     Felipe Galarce at INRIA/WIAS/PUCV

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

#include<wave_complex.hpp>

void WaveComplex::initialize(Parameters parameters, const Geometry & geometry, const Boundary & boundary){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "Elasticity frecuency: Initializing" << endl;
  par = parameters;
  geo = geometry;
  bd = boundary;
  m_verbose = par.verbose();

  nbDofVar.resize(par.nbVariables());
  nbDofs = 0;
  for (int i = 0; i < par.nbDofsPerNode().size(); i++){
    nbDofs += par.nbDofsPerNode()[i]*geo.nbVertices;
    nbDofVar[i] = par.nbDofsPerNode()[i]*geo.nbVertices;
  }
  nbVertices = geo.nbVertices;
  nbNodesPerElement = geo.dimension()+1; 

  mat(A, nbDofs, nbDofs); 
  code = MatGetOwnershipRange(A, &m, &n); CHKERR(code);

  fem.initialize(par, geo, bd, A);
}

void WaveComplex::finalize(){
  MatDestroy(&A);
  KSPDestroy(&ksp);
}

Mat WaveComplex::assembleLHS(){

  int m,n;
  code = MatGetOwnershipRange(A, &m, &n); CHKERR(code);
  if (m_world_rank == 0) cout << "WaveComplex: Assembling discretization matrix." << endl;
  for (int partId = 0; partId < geo.elements().size(); partId++){
    for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on elements */

      if (feId % (geo.elements()[partId].size() / m_verbose) == 0){ if (m_world_rank == 0) cout << "  Elementary matrix for element: " << feId << "/" << geo.elements()[partId].size() - 1 << endl; }

      vector<int> simplex = geo.elements()[partId][feId];

      fem.setSimplex(simplex);

      for (int var = 0; var < par.nbVariables(); var++){
        fem.u_dot_v_scalar(-par.frecuency()*par.frecuency(), var, var);

        for (int comp = 0; comp < geo.dimension(); comp++){
          fem.du_x_dv_y(par.wave_velocity()[0]*par.wave_velocity()[0], comp, comp, var, var);
        }
      }

      for (int var = 0; var < par.nbVariables()/2; var++){
        fem.u_dot_v_scalar(-par.frecuency()*par.damping(), var,                 geo.dimension()+var);
        fem.u_dot_v_scalar(-par.frecuency()*par.damping(), geo.dimension()+var, var                );
      }
    }
  }
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERR(code);

  /* LHS static matrix */
  double normLHS = norm(A);
  if (m_world_rank == 0) cout << "Poro complex: Norm LHS = " << normLHS << endl;
  
  return A;
}

Vec WaveComplex::assembleRHS(){
  Vec b = zeros(nbDofs);
  return b;
}

void WaveComplex::setSolver(){
  if (m_world_rank == 0) cout << "KSP: Configuring" << endl;
  configureKSP(ksp, par.solver(), par.preconditioner(), par.monitorKSP(), par.use_solution_as_guess_KSP(), par.reuse_preconditioner());
}

void WaveComplex::setLHS(Mat A, Mat P){
  if (m_world_rank == 0) cout << "KSP: Seting operators" << endl;
  if (P == NULL){
    code = KSPSetOperators(ksp, A, A); CHKERR(code); 
  } else {
    code = KSPSetOperators(ksp, A, P); CHKERR(code); 
  }
  code = KSPSetUp(ksp); CHKERR(code);
}

void WaveComplex::solve(Vec b, Vec u){
  if (m_world_rank == 0) cout << "WaveComplex: Solving linear system" << endl;
  code = KSPSolve(ksp, b, u); CHKERR(code); 
  int its;
  double rnorm;
  code = KSPGetResidualNorm(ksp, &rnorm); CHKERR(code);
  code = KSPGetIterationNumber(ksp, &its); CHKERR(code);
  if (par.solver() == "gmres"){
    if (m_world_rank == 0) cout << "WaveComplex: " << its << " iterations. Residual norm = " << rnorm << endl;
  }
}
