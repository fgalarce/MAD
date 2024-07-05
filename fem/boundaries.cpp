/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2024,
    
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

#include<boundaries.hpp>

void Boundary::initialize(Parameters parameters, const Geometry & geometry){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "Boundary: Initializing." << endl;
  par = parameters;
  geo = geometry;
  m_verbose = par.verbose();
  m_dimension = geo.dimension();
  nbDofs = 0;
  for (int i = 0; i < par.nbVariables(); i++){
    nbDofs = nbDofs + par.nbDofsPerNode()[i]*geo.nbVertices;
  }

  m_NeumannBC_weak = zeros(nbDofs);
  if (m_dimension > 1){
  
    feBD.initialize(par, m_dimension);
  
    /* mass matrix on the boundary */
    mat(m_M, nbDofs, nbDofs);
    m_NeumannMass.resize(geo.boundaryNodes().size());
    for (int i = 0; i < geo.boundaryNodes().size(); i++){
      mat(m_NeumannMass[i], nbDofs, nbDofs);
    }
  
    if (m_dimension == 3){
      computeCorners();
    }
    if (par.nbDofsPerNode()[0] > 1){
      assembleBoundaryNormals("contiguous");
      assembleBoundaryMass("contiguous");
    } else {
      assembleBoundaryMass("splitted", geo.dimension());
      assembleBoundaryNormals("splitted");
    }
  
    m_scalarMass.resize(geo.boundaryNodes().size()); 
 
    assembleBoundaryMassScalar(); // bifurcate this only if field is scalar
  
    if (par.hemodynamics()){
      if (m_world_rank == 0) cout << "- - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
      for (int i = 0; i < geo.elementsBD().size(); i++){
        m_paraboloid_center.push_back(geo.getCenter( geo.elementLabelsBD()[i]));
        m_paraboloid_radius.push_back(geo.getRadius( geo.elementLabelsBD()[i], m_paraboloid_center.back()));
        double size =  computeSize(geo.elementLabelsBD(i));
        if (m_world_rank == 0) cout << "Boundary: label = " << geo.elementLabelsBD()[i] << endl;
        if (m_world_rank == 0) cout << " center= " << m_paraboloid_center.back() << endl;
        if (m_world_rank == 0) cout << " radius = " << m_paraboloid_radius.back() << endl;
        if (m_world_rank == 0) cout << " size = " << size << endl;
      }
      if (par.inflow_data() != ""){
        m_inletU0 = importdata(par.inflow_data());
      }
      if (m_world_rank == 0) cout << "- - - - - - - - - - - - - - - - - - - - - - - - -" << endl;
    }
  
  }

}

void Boundary::block(Mat & A, Vec & rhs, string method){
  block(A, method);
  block(rhs, method);
}

void Boundary::block(Mat & A, string method){
  if (method == "brut"){
    for (int i=0;i<m_bc.size();i++){
      code = MatZeroRows(A, m_idDofsBC[i].size(), &m_idDofsBC[i][0], 1.0, NULL, NULL); CHKERR(code);
    }
  } else if (method == "symmetric"){
    int high, low;
    code = MatGetOwnershipRange(A, &low, &high); CHKERR(code);
    for (int i=0;i<m_bc.size();i++){
      for (int j = 0; j < m_idDofsBC[i].size(); j++){
        if (m_idDofsBC[i][j] >= low and m_idDofsBC[i][j] < high){
          MatSetValue(A, m_idDofsBC[i][j], m_idDofsBC[i][j], 1.0e+30, INSERT_VALUES);
        }
      }
    }
  }
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERR(code);
}

void Boundary::block(Vec & rhs, string method){
  /* Apply Neumann */
  code = VecAssemblyBegin(rhs); CHKERR(code);
  code = VecAssemblyEnd(rhs); CHKERR(code);
  code = VecAXPY(rhs, 1.0, m_NeumannBC_weak); CHKERR(code);
  /* Apply Dirichlet */
  double coef = 1.0;
  if (method == "symmetric"){
    coef = 1.0e+30;
  }
  int high, low;
  code = VecGetOwnershipRange(rhs, &low, &high); CHKERR(code);
  for (int i = 0; i < m_bc.size(); i++){
    for (int j = 0; j < m_idDofsBC[i].size(); j++){
      if (m_idDofsBC[i][j] >= low and m_idDofsBC[i][j] < high){
        code = VecSetValue(rhs, m_idDofsBC[i][j], m_bc[i][j]*coef, INSERT_VALUES); CHKERR(code);
      }
    }
  }
  code = VecAssemblyBegin(rhs); CHKERR(code);
  code = VecAssemblyEnd(rhs); CHKERR(code);
}

void Boundary::reset(){
  m_bc.clear();
  m_idDofsBC.clear();
  VecZeroEntries(m_NeumannBC_weak);
}

void Boundary::time(double t){
  if (m_time < t){
    m_iteration = m_iteration + 1;
  }
  m_time = t;
  reset();
}

void Boundary::Dirichlet(int bdLabel, vector<double> (*bc_func) (vector<double>,double,Parameters), int variableLabel){
  if (m_dimension > 1){
    Dirichlet_2D_or_3D(bdLabel, bc_func, variableLabel);
  } else {

    int offset = 0; /* to select variable */
    for (int i = 0; i < variableLabel; i++){
      offset = offset + par.nbDofsPerNode()[i]*geo.nbVertices;
    }
    vector<double> BC(1);
    vector<int> iBC(1);
    if (bdLabel == 0){
      BC[0] = bc_func(geo.coordinates().front(), m_time, par)[0];
      iBC[0] = offset;
    } else if (bdLabel == 1){
      BC[0] = bc_func(geo.coordinates().back(), m_time, par)[0];
      iBC[0] = offset + par.nbDofsPerNode()[variableLabel]*geo.nbVertices - 1;
    } else {
      if (m_world_rank == 0) cout << "ERROR: bd labels tolerated for one dimensional problems are either 0 (left) or 1 (right)." << endl; 
      exit(1);
    }
    m_bc.push_back(BC);
    m_idDofsBC.push_back(iBC);
    BC.clear();
    iBC.clear();
  }
}

void Boundary::Dirichlet_2D_or_3D(int bdLabel, vector<double> (*bc_func) (vector<double>,double,Parameters), int variableLabel){
  int nbDofsBC = par.nbDofsPerNode()[variableLabel]*geo.boundaryNodes(bdLabel).size();
  vector<double> BC(nbDofsBC);
  vector<int> iBC(nbDofsBC);
  int offset = 0; /* to select variable */
  for (int i = 0; i < variableLabel; i++){
    offset = offset + par.nbDofsPerNode()[i]*geo.nbVertices;
  }
  for (int i = 0; i < geo.boundaryNodes(bdLabel).size(); i++){
    vector<double> u_i = bc_func(geo.coordinates(geo.boundaryNodes(bdLabel)[i]), m_time, par);
    for (int comp = 0; comp < par.nbDofsPerNode()[variableLabel]; comp++){
      iBC[par.nbDofsPerNode()[variableLabel]*i+comp] = offset + par.nbDofsPerNode()[variableLabel]*geo.boundaryNodes(bdLabel)[i]+comp;
      BC[par.nbDofsPerNode()[variableLabel]*i+comp] = u_i[comp];
    }
  }
  m_bc.push_back(BC);
  m_idDofsBC.push_back(iBC);
  BC.clear();
  iBC.clear();
  if (m_world_rank == 0 && !shut_up) cout << "  Boundary: Dirichlet condition for variable '" << par.variableName()[variableLabel] <<"' over sub-domain: " << bdLabel << ".\n";
  if (m_world_rank == 0 && !shut_up) cout << "  Boundary: Number of degrees of freedom modified: " << nbDofsBC << ".\n";
}

void Boundary::DirichletParaboloid(int bdLabel, vector<double> (*bc_func) (vector<double>,double,Parameters), int variableLabel){
  int nbDofsBC = par.nbDofsPerNode()[variableLabel]*geo.boundaryNodes(bdLabel).size();
  vector<double> BC(nbDofsBC);
  vector<int> iBC(nbDofsBC);
  int offset = 0; /* to select variable */
  for (int i = 0; i < variableLabel; i++){
    offset = offset + par.nbDofsPerNode()[i]*geo.nbVertices;
  }
  for (int i = 0; i < geo.boundaryNodes(bdLabel).size(); i++){
    vector<double> point = geo.coordinates(geo.boundaryNodes(bdLabel)[i]);
    double paraboloid = (1.0 - dot(point-paraboloid_center(bdLabel), point-paraboloid_center(bdLabel)) / ( paraboloid_radius(bdLabel)*paraboloid_radius(bdLabel)));
    for (int comp = 0; comp < par.nbDofsPerNode()[variableLabel]; comp++){
      iBC[par.nbDofsPerNode()[variableLabel]*i + comp] = offset + par.nbDofsPerNode()[variableLabel]*geo.boundaryNodes(bdLabel)[i]+comp;
      BC[par.nbDofsPerNode()[variableLabel]*i + comp] = bc_func(point, m_time, par)[comp] * paraboloid;
    }
  }
  m_bc.push_back(BC);
  m_idDofsBC.push_back(iBC);
  if (m_world_rank == 0) cout << "  Boundary: Dirichlet condition over sub-domain: " << bdLabel << ".\n";
  if (m_world_rank == 0) cout << "  Boundary: Number of degrees of freedom modified: " << nbDofsBC << ".\n";
}

void Boundary::DirichletNormalParaboloid(int bdLabel, double (*bc_func) (vector<double>,double,Parameters), int variableLabel){
  int nbDofsPerNode = par.nbDofsPerNode()[variableLabel];
  int nbDofsBC = nbDofsPerNode*geo.boundaryNodes(bdLabel).size();
  vector<double> BC(nbDofsBC);
  vector<int> iBC(nbDofsBC);
  int offset = 0; /* to select variable */
  for (int i = 0; i < variableLabel; i++){
    offset = offset + par.nbDofsPerNode()[i]*geo.nbVertices;
  }
  for (int i = 0; i < geo.boundaryNodes(bdLabel).size(); i++){
    int pointId = geo.boundaryNodes(bdLabel)[i];
    vector<double> point = geo.coordinates(pointId);
    double paraboloid = (1.0 - dot(point-paraboloid_center(bdLabel), point-paraboloid_center(bdLabel)) / ( paraboloid_radius(bdLabel)*paraboloid_radius(bdLabel)));
    for (int comp = 0; comp < nbDofsPerNode; comp++){
      iBC[nbDofsPerNode*i + comp] = offset + nbDofsPerNode*pointId + comp;
      BC[nbDofsPerNode*i + comp] = m_stl_normals_sequential[nbDofsPerNode*pointId+comp] * bc_func(point, m_time, par) * paraboloid;
    }
  }
  m_bc.push_back(BC);
  m_idDofsBC.push_back(iBC);
  if (m_world_rank == 0) cout << "  Boundary: Dirichlet condition over sub-domain: " << bdLabel << ".\n";
  if (m_world_rank == 0) cout << "  Boundary: Number of degrees of freedom modified: " << nbDofsBC << ".\n";
}

void Boundary::DirichletNormalParaboloid(int bdLabel, double (*bc_func) (vector<double>,double,Parameters), vector<int> variableLabel){
  int nbDofsBC = variableLabel.size()*geo.boundaryNodes(bdLabel).size();
  vector<double> BC(nbDofsBC);
  vector<int> iBC(nbDofsBC);
  int offset = 0; /* to select variable */
  for (int i = 0; i < variableLabel[0]; i++){
    offset = offset + par.nbDofsPerNode()[i]*geo.nbVertices;
  }
  
  for (int i = 0; i < geo.boundaryNodes(bdLabel).size(); i++){
    int pointId = geo.boundaryNodes(bdLabel)[i];
    vector<double> point = geo.coordinates(pointId);
    double paraboloid = (1.0 - dot(point-paraboloid_center(bdLabel), point-paraboloid_center(bdLabel)) / ( paraboloid_radius(bdLabel)*paraboloid_radius(bdLabel)));
    for (int var = 0; var < variableLabel.size(); var++){
      iBC[i + var*geo.boundaryNodes(bdLabel).size()] = offset + var*geo.nbVertices + pointId;
      BC[i + var*geo.boundaryNodes(bdLabel).size()] = m_stl_normals_sequential[pointId+var*geo.nbVertices] * bc_func(point, m_time, par) * paraboloid;
    }
  }
  m_bc.push_back(BC);
  m_idDofsBC.push_back(iBC);
  if (m_world_rank == 0) cout << "  Boundary: Dirichlet condition over sub-domain: " << bdLabel << ".\n";
  if (m_world_rank == 0) cout << "  Boundary: Number of degrees of freedom modified: " << nbDofsBC << ".\n";
}

void Boundary::Dirichlet(int bdLabel, Vec u, int variableLabel){
  int nbDofsBC = par.nbDofsPerNode()[variableLabel]*geo.boundaryNodes(bdLabel).size();
  vector<double> BC(nbDofsBC);
  vector<int> iBC(nbDofsBC);
  int offset = 0; /* to select variable */
  for (int i = 0; i < variableLabel; i++){
    offset = offset + par.nbDofsPerNode()[i]*geo.nbVertices;
  }
  for (int i = 0; i < geo.boundaryNodes(bdLabel).size(); i++){
    for (int comp = 0; comp < par.nbDofsPerNode()[variableLabel]; comp++){
      iBC[par.nbDofsPerNode()[variableLabel]*i+comp] = offset + par.nbDofsPerNode()[variableLabel]*geo.boundaryNodes(bdLabel)[i]+comp;
    }
  }
  code = VecGetValues(getSequential(u), iBC.size(), &range(iBC.size())[0], &BC[0]); CHKERR(code);
  m_bc.push_back(BC);
  m_idDofsBC.push_back(iBC);
  if (m_world_rank == 0) cout << "  Boundary: Dirichlet condition for variable '" << par.variableName()[variableLabel] <<"' over sub-domain: " << bdLabel << ".\n";
  if (m_world_rank == 0) cout << "  Boundary: Number of degrees of freedom modified: " << nbDofsBC << ".\n";
}

void Boundary::Neumann(int bdLabel, vector<double> (*bc_func) (vector<double>, double, Parameters), int variableLabel){
  if (m_dimension > 1){
    Neumann_2D_or_3D(bdLabel, bc_func, variableLabel);
  } else {
    exit(1);
  }
}

void Boundary::Neumann_2D_or_3D(int bdLabel, vector<double> (*bc_func) (vector<double>, double, Parameters), int variableLabel){
  int nbDofsBC = par.nbDofsPerNode()[variableLabel]*geo.boundaryNodes(bdLabel).size();
  vector<double> BC(nbDofsBC);
  vector<int> iBC(nbDofsBC);
  int offset = 0; /* to select variable */
  for (int i = 0; i < variableLabel; i++){
    offset = offset + par.nbDofsPerNode()[i]*geo.nbVertices;
  }
  for (int i = 0; i < geo.boundaryNodes(bdLabel).size(); i++){
    vector<double> u_i = bc_func(geo.coordinates(geo.boundaryNodes(bdLabel)[i]), m_time, par);
    for (int comp = 0; comp < par.nbDofsPerNode()[variableLabel]; comp++){
      BC[par.nbDofsPerNode()[variableLabel]*i+comp] = u_i[comp];
      iBC[par.nbDofsPerNode()[variableLabel]*i+comp] = offset + par.nbDofsPerNode()[variableLabel]*geo.boundaryNodes(bdLabel)[i]+comp;
    }
  }
  Vec NeumannBC = zeros(nbDofs);
  Vec mass_times_NeumannBC = zeros(nbDofs);
  code = VecSetValues(NeumannBC, nbDofsBC, &iBC[0], &BC[0], INSERT_VALUES); CHKERR(code);
  code = VecAssemblyBegin(NeumannBC); CHKERR(code);
  code = VecAssemblyEnd(NeumannBC); CHKERR(code);
  code = MatMult(NeumannMass(bdLabel), NeumannBC, mass_times_NeumannBC); CHKERR(code);
  code = VecAXPY(m_NeumannBC_weak, 1.0, mass_times_NeumannBC); CHKERR(code);
  double norm_bc_weak;
  code = VecNorm(mass_times_NeumannBC, NORM_2, &norm_bc_weak); CHKERR(code);

  if (m_world_rank == 0) cout << "  Boundary: Neumann condition for variable '" << par.variableName()[variableLabel] <<"' over sub-domain: " << bdLabel << ".\n";
  if (m_world_rank == 0) cout << "  Boundary: Number of degrees of freedom modified: " << nbDofsBC << ".\n";
  if (m_world_rank == 0) cout << "  Boundary: Neumann BC rhs norm (concentred force): " << norm_bc_weak << ".\n";
}

void Boundary::NeumannNormal(int bdLabel, double (*bc_func) (vector<double>, double, Parameters), int variableLabel){
  int nbDofsPerNode = par.nbDofsPerNode()[variableLabel];
  int nbDofsBC = nbDofsPerNode*geo.boundaryNodes(bdLabel).size();
  vector<double> BC(nbDofsBC);
  vector<int> iBC(nbDofsBC);
  int offset = 0; /* to select variable */
  for (int i = 0; i < variableLabel; i++){
    offset = offset + par.nbDofsPerNode()[i]*geo.nbVertices;
  }
  for (int i = 0; i < geo.boundaryNodes(bdLabel).size(); i++){
    int pointId = geo.boundaryNodes(bdLabel)[i];
    vector<double> point = geo.coordinates(pointId);
    double u_i = bc_func(point, m_time, par);
    for (int comp = 0; comp < nbDofsPerNode; comp++){
      iBC[nbDofsPerNode*i + comp] = offset + nbDofsPerNode*geo.boundaryNodes(bdLabel)[i]+comp;
      BC[nbDofsPerNode*i + comp] = m_stl_normals_sequential[nbDofsPerNode*pointId+comp] * u_i;
    }
  }
  Vec NeumannBC = zeros(nbDofs);
  Vec mass_times_NeumannBC = zeros(nbDofs);
  code = VecSetValues(NeumannBC, nbDofsBC, &iBC[0], &BC[0], INSERT_VALUES); CHKERR(code);
  code = VecAssemblyBegin(NeumannBC); CHKERR(code);
  code = VecAssemblyEnd(NeumannBC); CHKERR(code);
  code = MatMult(NeumannMass(bdLabel), NeumannBC, mass_times_NeumannBC); CHKERR(code);
  code = VecAXPY(m_NeumannBC_weak, 1.0, mass_times_NeumannBC); CHKERR(code);
  double norm_bc_weak;
  code = VecNorm(mass_times_NeumannBC, NORM_2, &norm_bc_weak); CHKERR(code);

  if (m_world_rank == 0) cout << "  Boundary: Neumann condition for variable '" << par.variableName()[variableLabel] <<"' over sub-domain: " << bdLabel << ".\n";
  if (m_world_rank == 0) cout << "  Boundary: Number of degrees of freedom modified: " << nbDofsBC << ".\n";
  if (m_world_rank == 0) cout << "  Boundary: Neumann BC rhs norm (concentred force): " << norm_bc_weak << ".\n";
}

void Boundary::NeumannNormalParaboloid(int bdLabel, double (*bc_func) (vector<double>,double,Parameters), int variableLabel){
  vector<double> BC(nbDofs);
  for (int i = 0; i < geo.coordinates().size(); i++){
    double paraboloid = (1.0 - dot(geo.coordinates()[i] - paraboloid_center(bdLabel),  geo.coordinates()[i] - paraboloid_center(bdLabel)) / (paraboloid_radius(bdLabel)*paraboloid_radius(bdLabel)));
    for (int comp = 0; comp < par.nbDofsPerNode()[variableLabel]; comp++){
      BC[par.nbDofsPerNode()[variableLabel]*i+comp] 
          = m_stl_normals_sequential[par.nbDofsPerNode()[variableLabel]*i+comp] * bc_func(geo.coordinates()[i], m_time, par) * paraboloid;
    }
  }
  Vec NeumannBC = petsc(BC);
  code = MatMult(NeumannMass(bdLabel), NeumannBC, m_NeumannBC_weak); CHKERR(code);
  double norm_force;
  code = VecNorm(m_NeumannBC_weak, NORM_2, &norm_force); CHKERR(code);
  if (m_world_rank == 0) cout << "  Boundary: Neumann condition for variable '" << par.variableName()[variableLabel] <<"' over sub-domain: " << bdLabel << ".\n";
  if (m_world_rank == 0) cout << "  Boundary: Force norm: " << norm_force << endl;
}

void Boundary::NeumannNormalConstant(int bdLabel, double bc_value){
//  Vec NeumannBC = zeros(nbDofs);
//  Vec neumannBC_weak = zeros(nbDofs);
//  code = VecAXPY(NeumannBC, bc_value, m_normals); CHKERR(code);
//  code = MatMult(NeumannMass(bdLabel), NeumannBC, neumannBC_weak); CHKERR(code);
//  code = VecAXPY(m_NeumannBC_weak, 1.0, neumannBC_weak);
//  double norm_bc_weak;
//  code = VecNorm(neumannBC_weak, NORM_2, &norm_bc_weak); CHKERR(code);
//  if (m_world_rank == 0) cout << "  Boundary: Neumann boundary condition over sub-domain: " << bdLabel << ".\n";
//  if (m_world_rank == 0) cout << "  Boundary: Force norm: " << norm_bc_weak << endl;
//  VecDestroy(&neumannBC_weak);
//  VecDestroy(&NeumannBC);

  Vec neumannBC_weak = zeros(nbDofs);
  int nbNodesPerElement = geo.dimension();
  if (m_world_rank == 0) cout << "Boundary: Assembling Neumann RHS. \n";

  for (vector<int> simplex : geo.elementsBD(bdLabel)){ /* loop on boundary elements */
    /* get finite element coordinates */
    vector<vector<double>> coordinates(nbNodesPerElement);
    for (int i = 0; i < nbNodesPerElement; i++){
      coordinates[i] = geo.coordinates()[simplex[i]];
    }
    feBD.setCoordinates(coordinates);
    feBD.computeNormal();

    /* set values from local to global */
    for (int comp = 0; comp < geo.dimension(); comp++){
      for (int i = 0; i < nbNodesPerElement; i++){
        double integral = 0;
        for (int qp = 0; qp < par.nbQuadraturePointsBD(); qp++){
          integral = integral + feBD.phi(qp, i) * feBD.weights()[qp] ;
        }
        integral = integral * feBD.normal()[comp] * feBD.sqrtDetJTJ();
        code = VecSetValue(neumannBC_weak, simplex[i]+comp*geo.nbVertices, integral, ADD_VALUES); CHKERR(code);
      }
    }
  }

  code = VecAssemblyBegin(neumannBC_weak); CHKERR(code);
  code = VecAssemblyEnd(neumannBC_weak); CHKERR(code);
  code = VecAXPY(m_NeumannBC_weak, bc_value, neumannBC_weak);
  double norm_bc_weak;
  code = VecNorm(neumannBC_weak, NORM_2, &norm_bc_weak); CHKERR(code);
  if (m_world_rank == 0) cout << "  Boundary: Neumann boundary condition over sub-domain: " << bdLabel << ".\n";
  if (m_world_rank == 0) cout << "  Boundary: Force norm: " << norm_bc_weak << endl;
  VecDestroy(&neumannBC_weak);
}

void Boundary::backflow(Vec u){
  if (m_world_rank == 0) cout << "Boundary: Back-flow stabilization." << endl;
  Vec u_sequential = getSequential(u);
  vector<double> u_seq = stl(u_sequential);
  for (int idBD : par.outlets()){
    for (int i : geo.boundaryNodes(idBD)){
      double stabilize = 0.0;
      for (int comp = 0; comp < geo.dimension(); comp++){
        stabilize = stabilize 
                  + u_seq[i + geo.nbNodes()*comp]
                  * m_stl_normals_sequential[i + geo.nbNodes()*comp];
      } 
      if (stabilize < 0){
        for (int comp = 0; comp < geo.dimension(); comp++){
          vecSetInsert(u, i + geo.nbNodes()*comp, 0.0);
        }
      }
    }
  } 
  code = VecAssemblyBegin(u); CHKERR(code);
  code = VecAssemblyEnd(u); CHKERR(code);
  VecDestroy(&u_sequential);
  u_seq.clear();
}

void Boundary::backflowNeumann(Vec u, int bdLabel){
  Vec stab_force = zeros(nbDofs);
  Vec u_sequential = getSequential(u);
  vector<double> u_seq = stl(u_sequential);
  for (int i : geo.boundaryNodes(bdLabel)){
    double normal_u = 0.0;
    for (int var = 0; var < par.nbVariables()-1; var++){
      normal_u += u_seq[i + geo.nbNodes()*var]
                * m_stl_normals_sequential[i + geo.nbNodes()*var];
    }
    for (int var = 0; var < par.nbVariables()-1; var++){
      vecSetInsert(stab_force,
                    i + geo.nbNodes()*var,
                    u_seq[i + geo.nbNodes()*var]*0.5*(normal_u - std::abs(normal_u)));
    }
  }

  code = VecAssemblyBegin(stab_force); CHKERR(code);
  code = VecAssemblyEnd(stab_force); CHKERR(code);
  VecScale(stab_force, par.density()*par.backflowCoeff());

  Vec NeumannBC_weak = zeros(nbDofs);
  code = MatMult(NeumannMass(bdLabel), stab_force, NeumannBC_weak); CHKERR(code);
  code = VecAXPY(m_NeumannBC_weak, 1.0, NeumannBC_weak); CHKERR(code);
  double norm_bc_weak;
  code = VecNorm(NeumannBC_weak, NORM_2, &norm_bc_weak); CHKERR(code);

  if (m_world_rank == 0) cout << "  Boundary: Neumann back-flow stabilization." << endl;
  if (m_world_rank == 0) cout << "  Boundary: Neumann boundary condition over sub-domain: " << bdLabel << ".\n";
  if (m_world_rank == 0) cout << "  Boundary: Force norm: " << norm_bc_weak << endl;

  VecDestroy(&u_sequential);
  VecDestroy(&stab_force);
  VecDestroy(&NeumannBC_weak);
  u_seq.clear();
}


void Boundary::backflowStrong(Vec u, int idBD){
  Vec u_sequential = getSequential(u);
  vector<double> u_seq = stl(u_sequential);
  for (int i : geo.boundaryNodes(idBD)){
    double stabilize = 0.0;
    for (int comp = 0; comp < geo.dimension(); comp++){
      stabilize = stabilize + u_seq[par.nbDofsPerNode()[0]*i+comp]*m_stl_normals_sequential[par.nbDofsPerNode()[0]*i+comp];
    } 
    if (stabilize < 0){
      for (int j : geo.boundaryNodes(idBD)){
        for (int comp = 0; comp < geo.dimension(); comp++){
          vecSetInsert(u, par.nbDofsPerNode()[0]*j+comp, 0.0);
        }
      }
      code = VecAssemblyBegin(u); CHKERR(code);
      code = VecAssemblyEnd(u); CHKERR(code);
      VecDestroy(&u_sequential);
      u_seq.clear();
      return;
    }
  }
}

void Boundary::backflowStrong(Vec u){
  if (m_world_rank == 0) cout << "Boundary: Back-flow stabilization." << endl;
  for (int idBD : par.outlets()){
    backflowStrong(u, idBD);
  } 
}

void Boundary::DirichletComp(int compBC, int bdLabel, vector<double> (*bc_func) (vector<double>,double,Parameters), int variableLabel){
  int nbDofsBC = geo.boundaryNodes(bdLabel).size();
  vector<double> BC(nbDofsBC);
  vector<int> iBC(nbDofsBC);
  int offset = 0; /* to select variable */
  for (int i = 0; i < variableLabel; i++){
    offset = offset + par.nbDofsPerNode()[i]*geo.nbVertices;
  }
  for (int i = 0; i < geo.boundaryNodes(bdLabel).size(); i++){
    vector<double> u_i = bc_func(geo.coordinates(geo.boundaryNodes(bdLabel)[i]), m_time, par);
    iBC[i] = offset + par.nbDofsPerNode()[variableLabel]*geo.boundaryNodes(bdLabel)[i]+compBC;
    BC[i] = u_i[compBC];
  }
  m_bc.push_back(BC);
  m_idDofsBC.push_back(iBC);
  if (m_world_rank == 0) cout << "  Boundary: " << compBC << "-component Dirichlet condition for variable '" << par.variableName()[variableLabel] <<"' over sub-domain: " << bdLabel << ".\n";
  if (m_world_rank == 0) cout << "  Boundary: Number of degrees of freedom modified: " << nbDofsBC << ".\n";
}

void Boundary::assembleBoundaryMass(string mode, int nbVar){

  code = MatAssemblyBegin(m_M, MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyEnd(m_M, MAT_FINAL_ASSEMBLY); CHKERR(code);
  int idElement = 0;
  int nbNodesPerElement = geo.dimension();
  if (m_world_rank == 0) cout << "Boundary: Assembling mass matrix on the boundary. \n";
  for (int partId = 0; partId < geo.elementsBD().size(); partId++){ /* loop on parts */
    for (vector<int> simplex : geo.elementsBD()[partId]){ /* loop on boundary elements */
      /* get finite element coordinates */
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int i = 0; i < nbNodesPerElement; i++){
        coordinates[i] = geo.coordinates()[simplex[i]];
      }
      feBD.setCoordinates(coordinates);

      /* set values from local to global */
      if (mode == "contiguous"){
        for (int i = 0; i < nbNodesPerElement; i++){
          for (int j = 0; j < nbNodesPerElement; j++){
            double Mij = feBD.mass(i,j);
            for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
              matSet(m_NeumannMass[partId], par.nbDofsPerNode()[0]*simplex[i]+comp, par.nbDofsPerNode()[0]*simplex[j]+comp, Mij);
            }
          }
        }
      } else if (mode == "splitted") {
        for (int i = 0; i < nbNodesPerElement; i++){
          for (int j = 0; j < nbNodesPerElement; j++){
            double Mij = feBD.mass(i,j);
            for (int var = 0; var < nbVar; var++){
              matSet(m_NeumannMass[partId], simplex[i]+var*geo.nbVertices, simplex[j]+var*geo.nbVertices, Mij);
            }
          }
        }
      }
      idElement++;
    }
    code = MatAssemblyBegin(m_NeumannMass[partId], MAT_FINAL_ASSEMBLY); CHKERR(code);
    code = MatAssemblyEnd(m_NeumannMass[partId], MAT_FINAL_ASSEMBLY); CHKERR(code);
    code = MatAXPY(m_M, 1.0, m_NeumannMass[partId], DIFFERENT_NONZERO_PATTERN); CHKERR(code);
  }
  code = MatAssemblyBegin(m_M, MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyEnd(m_M, MAT_FINAL_ASSEMBLY); CHKERR(code);

}

void Boundary::assembleBoundaryNormals(string mode){

  m_normals_P0.resize(geo.dimension()*geo.nbElementsBoundary());
  m_normals_P0comp.resize(geo.dimension());
  for (int comp = 0; comp < geo.dimension(); comp++){
    m_normals_P0comp[comp].resize(geo.nbElementsBoundary());
  }

  int idElement = 0;
  int nbNodesPerElement = geo.dimension();
  if (m_world_rank == 0) cout << "Boundary: Computing normals. \n";
  for (int partId = 0; partId < geo.elementsBD().size(); partId++){ /* loop on parts */
    if (m_world_rank == 0) cout << "  Part: " << geo.elementLabelsBD(partId) << endl;
    for (vector<int> simplex : geo.elementsBD()[partId]){ /* loop on boundary elements */
      /* get finite element coordinates */
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int i = 0; i < nbNodesPerElement; i++){
        coordinates[i] = geo.coordinates()[simplex[i]];
      }
      feBD.setCoordinates(coordinates);
      feBD.computeNormal();
      if (mode == "contiguous"){
        for (int comp = 0; comp < geo.dimension(); comp++){
          m_normals_P0[geo.dimension()*idElement + comp] = feBD.normal()[comp];
        }
      } else if (mode == "splitted") {
        for (int comp = 0; comp < geo.dimension(); comp++){
          m_normals_P0comp[comp][idElement] = feBD.normal()[comp];
        }
      }

      idElement++;
    }
  }

  INT interpolator;
  interpolator.initialize(par, geo);
  Vec normals;
  vector<Vec> normalsComp(geo.dimension());
  m_normals = vec(nbDofs);
  m_normals_only = vec(m_dimension*geo.nbVertices);
  code = VecZeroEntries(m_normals); CHKERR(code);
  if (mode == "contiguous"){
    normals = interpolator.interpolateP0_P1_boundary(m_normals_P0, "contiguous");
    vector<double> stl_normals_sequential = stl(normals);
    for (int i = 0; i < stl_normals_sequential.size(); i++){
      vecSet(m_normals, i, stl_normals_sequential[i]);
      vecSet(m_normals_only, i, stl_normals_sequential[i]);
    }
  } else if (mode == "splitted") {
    for (int comp = 0; comp < geo.dimension(); comp++){
      normalsComp[comp] = interpolator.interpolateP0_P1_boundary(m_normals_P0comp[comp]);
      vector<double> stl_normals_sequential = stl(normalsComp[comp]);
      for (int i = 0; i < stl_normals_sequential.size(); i++){
        vecSet(m_normals, i + comp*geo.nbVertices, stl_normals_sequential[i]);
        vecSet(m_normals_only, i + comp*geo.nbVertices, stl_normals_sequential[i]);
      }
    }
  }
  code = VecAssemblyBegin(m_normals); CHKERR(code);
  code = VecAssemblyEnd(m_normals); CHKERR(code);

  m_normals_sequential = getSequential(m_normals);
  m_stl_normals_sequential = stl(m_normals_sequential);

}

void Boundary::assembleBoundaryMassScalar(){
  int nbNodesPerElement = geo.dimension();
  if (m_world_rank == 0) cout << "Boundary: Assembling scalar mass matrix on the boundary. \n";
  for (int partId = 0; partId < geo.elementsBD().size(); partId++){ /* loop on parts */
    mat(m_scalarMass[partId], geo.nbVertices, geo.nbVertices);
    for (vector<int> simplex : geo.elementsBD()[partId]){ /* loop on boundary elements */
      /* get finite element coordinates */
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int i = 0; i < nbNodesPerElement; i++){
        coordinates[i] = geo.coordinates()[simplex[i]];
      }
      feBD.setCoordinates(coordinates);
      /* set values from local to global */
      for (int i = 0; i < nbNodesPerElement; i++){
        for (int j = 0; j < nbNodesPerElement; j++){
          double Mij = feBD.mass(i,j);
          matSet(m_scalarMass[partId], simplex[i], simplex[j], Mij);
        }
      }
    }
    code = MatAssemblyBegin(m_scalarMass[partId], MAT_FINAL_ASSEMBLY); CHKERR(code);
    code = MatAssemblyEnd(m_scalarMass[partId], MAT_FINAL_ASSEMBLY); CHKERR(code);
  }
}

double Boundary::flow(Vec u, int labelBD){
  double Q_bd;
  Vec m_times_v = vec(nbDofs);
  code = MatMult(NeumannMass(labelBD), m_normals, m_times_v); CHKERR(code);
  code = VecDot(u, m_times_v, &Q_bd); CHKERR(code);
  Q_bd = Q_bd * par.density();
  VecDestroy(&m_times_v);
  return Q_bd;
}

double Boundary::integral(Vec p, int labelBD){
  double int_bd;
  Vec m_times_v = vec(geo.nbVertices);
  Vec uno = vec(geo.nbVertices);
  code = VecSet(uno, 1.0); CHKERR(code);
  code = VecAssemblyBegin(uno); CHKERR(code);
  code = VecAssemblyEnd(uno); CHKERR(code);
  code = MatMult(scalarMass(labelBD), uno, m_times_v); CHKERR(code);
  code = VecDot(p, m_times_v, &int_bd); CHKERR(code);
  return int_bd;
}

double Boundary::computeSize(int bdLabel){
  Vec uno = zeros(geo.nbVertices);
  for (int i : geo.boundaryNodes(bdLabel)){
    code = VecSetValue(uno, i, 1.0, INSERT_VALUES); CHKERR(code);
  }
  code = VecAssemblyBegin(uno); CHKERR(code);
  code = VecAssemblyEnd(uno); CHKERR(code);
  return integral(uno, bdLabel);
}

void Boundary::reflection(Vec u, int bdLabel){
  if (m_dimension == 1){
    assert(bdLabel == 0 or bdLabel == 1);
    int bdId;
    if (bdLabel == 0) {
      bdId = 0;
    } else {
      bdId = geo.nbVertices-1;
    }
    double bdValue = -stl(u)[bdId]; 

    vector<double> BC(1, bdValue);
    vector<int> iBC(1, bdId);
    m_bc.push_back(BC);
    m_idDofsBC.push_back(iBC);
    BC.clear();
    iBC.clear();
  } else {
    if (m_world_rank == 0) cout << "Boundary: reflexive bc not yet implemented for dimension greater than one." << endl;
    exit(1);
  }
}

void Boundary::computeCorners(){

  int nbNodesPerElement = geo.dimension();
  if (m_world_rank == 0) cout << "Boundary: Computing corners. \n";
  for (int partId = 0; partId < geo.boundaryNodes().size(); partId++){
    for (int i : geo.boundaryNodes()[partId]){
      vector<int> partIds = geo.getPointStarBD(i)[1]; /* If all elements have the same partId it means the point is not a corner*/
      for (int j = 0; j < partIds.size()-1; j++){
        if (partIds[j] != partIds[j+1]){
          m_corners.push_back(i);
          break;
        }
      }
    }
  }
}
