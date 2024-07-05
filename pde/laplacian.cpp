/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2023,
    
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

#include<laplacian.hpp>

void Laplacian::initialize(Parameters parameters, const Geometry & geometry, const Boundary & boundary){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "Laplacian: Initializing" << endl;
  par = parameters;
  geo = geometry;
  bd = boundary;
  fe.initialize(par, geo.dimension());
  m_verbose = par.verbose();
  nbDofs = 0;

  nbVertices = geo.nbVertices;
  nbDofs = nbVertices * par.nbVariables();

  nbNodesPerElement = geo.dimension()+1; 
  mat(A, nbDofs, nbDofs); 
  mat(M, nbDofs, nbDofs); 
  code = MatGetOwnershipRange(A, &m, &n); CHKERR(code);

  femStat.initialize(par, geo, bd, A);
  femMass.initialize(par, geo, bd, M);
}

void Laplacian::finalize(){
  KSPDestroy(&ksp);
}

Mat Laplacian::assembleLHS(){
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
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERR(code);
  double normLHS = norm(A);
  if (m_world_rank == 0) cout << "Laplacian: Norm LHS = " << normLHS << endl;
  return A;
}

void Laplacian::setSolver(){
  if (m_world_rank == 0) cout << "KSP: Configuring" << endl;
  configureKSP(ksp, par);
}

void Laplacian::setLHS(Mat C){
  if (m_world_rank == 0) cout << "KSP: Seting operators" << endl;
  code = KSPSetOperators(ksp, C, C); CHKERR(code); 
  code = KSPSetUp(ksp); CHKERR(code);
}

void Laplacian::solve(Vec b, Vec u){
    if (m_world_rank == 0) cout << "Laplacian: Solving linear system" << endl;
    code = KSPSolve(ksp, b, u); CHKERR(code); 
    int its;
    double rnorm;
    code = KSPGetResidualNorm(ksp, &rnorm); CHKERR(code);
    code = KSPGetIterationNumber(ksp, &its); CHKERR(code);
    if (par.solver() == "gmres"){
      if (m_world_rank == 0) cout << "Laplacian: " << its << " iterations. Residual norm = " << rnorm << endl;
    }
}
