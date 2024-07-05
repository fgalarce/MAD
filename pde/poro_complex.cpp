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

#include<poro_complex.hpp>

void PoroComplex::initialize(Parameters parameters, const Geometry & geometry, const Boundary & boundary, const IO & inputOutput){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "Poro-elasticity frecuency: Initializing" << endl;
  par = parameters;
  geo = geometry;
  bd = boundary;
  io = inputOutput;
  m_verbose = par.verbose();

  assert(par.nbDofsPerNode().size() == 6);
  nbDofVar.resize(6);
  nbDofs = 0;
  for (int i = 0; i < par.nbDofsPerNode().size(); i++){
    nbDofs += par.nbDofsPerNode()[i]*geo.nbVertices;
    nbDofVar[i] = par.nbDofsPerNode()[i]*geo.nbVertices;
  }
  nbVertices = geo.nbVertices;
  nbNodesPerElement = geo.dimension()+1; 

  m_uAnal = zeros(nbDofs);
  mat(A, nbDofs, nbDofs); 
  code = MatGetOwnershipRange(A, &m, &n); CHKERR(code);

  calculus.initialize(par, geo, bd);
  fem.initialize(par, geo, bd, A);
}

void PoroComplex::finalize(){
  MatDestroy(&A);
  KSPDestroy(&ksp);
}

Mat PoroComplex::assembleLHS(){

  double Skemptom = 0.99;
  double BiotWillis = 1.0;

  int m,n;
  code = MatGetOwnershipRange(A, &m, &n); CHKERR(code);
  if (m_world_rank == 0) cout << "PEDup: Assembling discretization matrix." << endl;
  for (int partId = 0; partId < geo.elements().size(); partId++){

    double lambda = par.youngModulus()[partId] * par.poissonRatio()[partId] / ((1.0 + par.poissonRatio()[partId])*(1.0 - 2.0*par.poissonRatio()[partId]));
    double mu = par.youngModulus()[partId]/(2.0*(1.0+par.poissonRatio()[partId]));
    double bulk  = par.youngModulus()[partId] / (3.0* ( 1.0 - 2.0*par.poissonRatio()[partId]));
    double permeability = par.permeability()[partId];

    if (m_world_rank == 0) cout << "  Assembling elementary matrices, part " << partId << " (label " << geo.elementLabels(partId) << ")" << endl;
    if (m_world_rank == 0) cout << "  mu = " << mu << endl;
    if (m_world_rank == 0) cout << "  lambda = " << lambda << endl;
    if (m_world_rank == 0) cout << "  Bulk = " << bulk << endl;
    if (m_world_rank == 0) cout << "  Permeability = " << permeability << endl;

    for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on elements */

      if (feId % (geo.elements()[partId].size() / m_verbose) == 0){ if (m_world_rank == 0) cout << "  Elementary matrix for element: " << feId << "/" << geo.elements()[partId].size() - 1 << endl; }

      vector<int> simplex = geo.elements()[partId][feId];

      fem.setSimplex(simplex);

      /* u : v  */
      fem.u_dot_v(-1.0 * par.frecuency()*par.frecuency() * par.density(), 0, 0);
      fem.u_dot_v(-1.0 * par.frecuency()*par.frecuency() * par.density(), 1, 1);

      /* \epsilon(u) : \epsilon(v)  */
      fem.epsilon_u_epsilon_v(2.0*mu, 0, 0);
      fem.epsilon_u_epsilon_v(2.0*mu, 1, 1);

      /* lambda* p \cdot r */
      fem.u_dot_v( (1.0/lambda + Skemptom/BiotWillis), 3, 2);
      fem.u_dot_v(-(1.0/lambda + Skemptom/BiotWillis), 2, 3);

      /* \grad (p) : \grad (q)  */
      fem.grad_u_grad_v(permeability/par.viscosity()/par.frecuency()/BiotWillis, 2, 2);
      fem.grad_u_grad_v(permeability/par.viscosity()/par.frecuency()/BiotWillis, 3, 3);

      /* -phi \cdot div v */
      fem.u_div_v(-1.0, 0, 4);
      fem.u_div_v(-1.0, 1, 5);

      /* xi \cdot div u_ */
      fem.div_u_v(1.0, 4, 0);
      fem.div_u_v(1.0, 5, 1);

      /* phi : xi */
      fem.u_dot_v(1.0/lambda, 4, 4);
      fem.u_dot_v(1.0/lambda, 5, 5);

      /* -1/lambda phi \cdot q */
      fem.u_dot_v(-1.0/lambda, 3, 4);
      fem.u_dot_v(1.0/lambda, 2, 5);

      /* -1/lambda p \cdot xi */
      fem.u_dot_v(-1.0/lambda, 4, 2);
      fem.u_dot_v(-1.0/lambda, 5, 3);

      // ---- stabilizer ----  //
      double tau = fem.feSize()*fem.feSize();
      double delta_1 = 1/par.frecuency()/par.frecuency();
      double delta_2 = 1.0;
      /* \omega^4 \rho^2 (u,v) */
      fem.u_dot_v(delta_1*tau*par.density()*par.density()*pow(par.frecuency(),4), 0, 0);
      fem.u_dot_v(delta_1*tau*par.density()*par.density()*pow(par.frecuency(),4), 1, 1);

      /* -4 \omega^2 \rho \mu (\nabla u, \nabla v) */
      fem.grad_u_grad_v(delta_1*tau*-4.0*par.frecuency()*par.frecuency()*par.density()*mu, 0, 0);
      fem.grad_u_grad_v(delta_1*tau*-4.0*par.frecuency()*par.frecuency()*par.density()*mu, 1, 1);

      /* 2 \omega^2 \rho (div(u), div(v)) */
      fem.div_u_div_v(delta_1 * tau*2*par.frecuency()*par.frecuency()*par.density(), 0, 0);
      fem.div_u_div_v(delta_1 * tau*2*par.frecuency()*par.frecuency()*par.density(), 0, 0);

      /* \omega^2 * \rho (div(u), xi) */
      fem.div_u_v(delta_1*tau * par.frecuency() * par.frecuency() * par.density(), 4, 0); 
      fem.div_u_v(delta_1*tau * par.frecuency() * par.frecuency() * par.density(), 5, 1); 

      /* \omega^2 * \rho (phi, div(v))  */
      fem.u_div_v(delta_1 * tau * par.frecuency() * par.frecuency() * par.density(), 0, 4);
      fem.u_div_v(delta_1 * tau * par.frecuency() * par.frecuency() * par.density(), 1, 5);

      /* \nabla phi, \nabla xi */
      fem.grad_u_grad_v(delta_2 * 1.0/BiotWillis * 1.0/par.viscosity() * 1.0/par.frecuency() * tau, 4, 4);
      fem.grad_u_grad_v(delta_2 * 1.0/BiotWillis * 1.0/par.viscosity() * 1.0/par.frecuency() * tau, 5, 5);

    }
  }
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERR(code);

  /* LHS static matrix */
  double normLHS = norm(A);
  if (m_world_rank == 0) cout << "Poro complex: Norm LHS = " << normLHS << endl;
  
  return A;
}

Vec PoroComplex::assembleRHS(vector<double> (*manufactured_solution) (vector<double>,double,Parameters), int idVariable){
  Vec b = zeros(nbDofs);
  int offset = 0;
  for (int i = 0; i < idVariable; i++){
    offset += par.nbDofsPerNode()[i]*geo.nbVertices;
  }
  /* assemble manufactured solution */
  Vec manuU = vec(nbDofs);
  manufacture(manuU, offset, manufactured_solution, idVariable);
  code = VecAssemblyBegin(manuU); CHKERR(code);
  code = VecAssemblyEnd(manuU); CHKERR(code);
  MatMult(A, manuU, b);
  code = VecAXPY(m_uAnal, 1.0, manuU); CHKERR(code);
  return b;
}

void PoroComplex::manufacture(Vec sol, int offset, vector<double> (*analytic) (vector<double>,double,Parameters), int idVariable){
  for (int i = 0; i < geo.coordinates().size(); i++){
    vector<double> xx = geo.coordinates()[i];
    vector<double> U = analytic(xx, 0.0, par);
    for (int j = 0; j < par.nbDofsPerNode()[idVariable]; j++){
      double iU = par.nbDofsPerNode()[idVariable]*i+j;
      vecSet(sol, offset + iU, U[j]);
    }
  }
}

Vec PoroComplex::assembleRHS(){
  Vec b = zeros(nbDofs);
  return b;
}

void PoroComplex::setSolver(){
  if (m_world_rank == 0) cout << "KSP: Configuring" << endl;
//  configureKSP(ksp, par.solver(), par.preconditioner(), par.monitorKSP(), par.use_solution_as_guess_KSP(), par.reuse_preconditioner());
  configureKSP(ksp, par);
}

void PoroComplex::setLHS(Mat A, Mat P){
  if (m_world_rank == 0) cout << "KSP: Seting operators" << endl;
  if (P == NULL){
    code = KSPSetOperators(ksp, A, A); CHKERR(code); 
  } else {
    code = KSPSetOperators(ksp, A, P); CHKERR(code); 
  }
  code = KSPSetUp(ksp); CHKERR(code);
}

void PoroComplex::solve(Vec b, Vec u){
  if (m_world_rank == 0) cout << "PoroComplex: Solving linear system" << endl;
  code = KSPSolve(ksp, b, u); CHKERR(code); 
  int its;
  double rnorm;
  code = KSPGetResidualNorm(ksp, &rnorm); CHKERR(code);
  code = KSPGetIterationNumber(ksp, &its); CHKERR(code);
  if (par.solver() == "gmres"){
    if (m_world_rank == 0) cout << "PoroComplex: " << its << " iterations. Residual norm = " << rnorm << endl;
  }
}

/* Assemble U and P ground truths. Return manufactured RHS with phi_analytic */
Vec PoroComplex::computePseudoPressure(vector<double> (*Ure) (vector<double>,double,Parameters),
                                       vector<double> (*Uim) (vector<double>,double,Parameters), 
                                       vector<double> (*Pre) (vector<double>,double,Parameters),
                                       vector<double> (*Pim) (vector<double>,double,Parameters)){

  /* assemble div u */
  Vec u_re = vec(nbDofVar[0]);
  Vec u_im = vec(nbDofVar[1]);
  
  manufacture(u_re, 0, Ure, 0);
  manufacture(u_im, 0, Uim, 1);

  code = VecAssemblyBegin(u_re); CHKERR(code);
  code = VecAssemblyEnd(u_re); CHKERR(code);
  code = VecAssemblyBegin(u_im); CHKERR(code);
  code = VecAssemblyEnd(u_im); CHKERR(code);

  Vec divuRe = calculus.divergence(u_re);
  Vec divuIm = calculus.divergence(u_im);

  double lambda = par.youngModulus()[0] * par.poissonRatio()[0] / ((1.0 + par.poissonRatio()[0])*(1.0 - 2.0*par.poissonRatio()[0]));
  /* p - lambda div u */
  Vec phi_re = vec(nbDofVar[4]);
  Vec phi_im = vec(nbDofVar[5]);
  manufacture(phi_re, 0, Pre, 2);
  manufacture(phi_im, 0, Pim, 3);
  code = VecAssemblyBegin(phi_re); CHKERR(code);
  code = VecAssemblyEnd(phi_re); CHKERR(code);
  code = VecAssemblyBegin(phi_im); CHKERR(code);
  code = VecAssemblyEnd(phi_im); CHKERR(code);
  code = VecAXPY(phi_re, -lambda, divuRe); CHKERR(code);
  code = VecAXPY(phi_im, -lambda, divuIm); CHKERR(code);

  Vec phi = vec(nbDofs);
  vector<double> phiRe_seq = stl(getSequential(phi_re));
  vector<double> phiIm_seq = stl(getSequential(phi_im));
  for (int i = 0; i < geo.coordinates().size(); i++){
    int iPhi_re = geo.nbVertices*par.nbDofsPerNode()[0]*2 + geo.nbVertices*2 + i;
    int iPhi_im = geo.nbVertices*par.nbDofsPerNode()[0]*2 + geo.nbVertices*3 + i;
    VecSetValue(phi, iPhi_re, phiRe_seq[i], INSERT_VALUES);
    VecSetValue(phi, iPhi_im, phiIm_seq[i], INSERT_VALUES);
  }
  code = VecAssemblyBegin(phi); CHKERR(code);
  code = VecAssemblyEnd(phi); CHKERR(code);

  Vec b = vec(nbDofs);
  MatMult(A, phi, b);
  code = VecAXPY(m_uAnal, 1.0, phi); CHKERR(code);
  return b;
}

double PoroComplex::innerProduct(Vec u, Vec v, Geometry & geometry){
  int u_size; double ip;
  code = VecGetSize(u, &u_size); CHKERR(code);
  Vec Mv = vec(u_size);
  if (u_size == par.nbDofsPerNode()[0]*geometry.nbVertices){
    code = MatMult(m_ip[0], v, Mv); CHKERR(code);
  } else if (u_size == par.nbDofsPerNode()[2]*geometry.nbVertices){
    code = MatMult(m_ip[1], v, Mv); CHKERR(code);
  } else {
    errorMessage("innerProduct", "Wrong vector size");
  }
  code = VecDot(Mv, u, &ip); CHKERR(code);
  code = VecDestroy(&Mv); CHKERR(code);

  return ip;
}

void PoroComplex::assembleInnerProduct(Geometry & geometry){

  m_ip.resize(2);
  int nbDofsScalar = par.nbDofsPerNode()[2]*geometry.nbVertices;
  int nbDofsVector = par.nbDofsPerNode()[0]*geometry.nbVertices;
  mat(m_ip[0], nbDofsVector, nbDofsVector);
  mat(m_ip[1], nbDofsScalar, nbDofsScalar);
  fe.initialize(par, geometry.dimension());

  for (int partId = 0; partId < geometry.elements().size(); partId++){
    if (m_world_rank == 0) cout << "Navier Stokes: Assembling inner product matrix. Part " << partId << endl;
    for (int feId = 0; feId < geometry.elements()[partId].size(); feId++){ /* loop on elements */

      vector<int> simplex = geometry.elements()[partId][feId];
      /* get finite element coordinates */
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int i = 0; i < nbNodesPerElement; i++){
        coordinates[i] = geometry.coordinates()[simplex[i]];
      }

      fe.setCoordinates(coordinates);
      fe.computeSize();

      if (feId % (geometry.elements()[partId].size() / m_verbose) == 0){
        if (m_world_rank == 0) cout << "  Elementary matrix for element: " << feId << "/" << geometry.elements()[partId].size() - 1 << endl;
      }
      /* Assemble elementary matrices */
      for (int j = 0; j < nbNodesPerElement; j++){
        for (int i = 0; i < nbNodesPerElement; i++){
          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
            if (par.nbDofsPerNode()[0]*simplex[i]+comp >= m && par.nbDofsPerNode()[0]*simplex[i]+comp < n){
              double Kij = fe.stiffness(i,j);
              double Mij = fe.mass(i,j);
              code = MatSetValue(m_ip[0], par.nbDofsPerNode()[0]*simplex[i]+comp, par.nbDofsPerNode()[0]*simplex[j]+comp, Kij+Mij, ADD_VALUES); CHKERR(code);
            }
          }
          if (simplex[i] >= m && simplex[i] < n){
            double Mij = fe.mass(i,j);
            code = MatSetValue(m_ip[1], simplex[i], simplex[j], Mij, ADD_VALUES); CHKERR(code); /* B^T */ 
          }
        }
      }
    }
  }
  code = MatAssemblyBegin(m_ip[0], MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(m_ip[0], MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyBegin(m_ip[1], MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(m_ip[1], MAT_FINAL_ASSEMBLY); CHKERR(code);
}
