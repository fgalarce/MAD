/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2020,
    
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

#include<calculus.hpp>

void Calculus::initialize(Parameters parameters, const Geometry & geometry, const Boundary & boundary){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "CALCULUS: Initializing." << endl;
  par = parameters;
  geo = geometry;
  bd = boundary;
  m_dimension = geo.dimension();
  m_interpolator.initialize(par, geo);
  assert(m_dimension == 1 or m_dimension == 2 or m_dimension == 3);
//  m_nbDofs = par.nbDofsPerNode()[0]*io.nbVertices();
  m_nbVertices = geo.nbVertices;
  fe.initialize(par, geo.dimension());
  feBD.initialize(par, geo.dimension());
}

Vec Calculus::gradient(Vec u, int partId){
  if (m_world_rank == 0) cout << "CALCULUS: Computing gradient." << endl;
  int vecSize;
  VecGetSize(u, &vecSize);
  if (vecSize == geo.nbVertices){
    return gradientScalar(u, partId);
  } else {
    errorMessage("Calculus::gradient", "Can only compute gradients for vector with size equal to the number of nodes in the mesh");  
  }
}

Vec Calculus::gradientScalar(Vec u, int bdLabel){

  vector<vector<double>> G = transpose(fe.grad_phi());
  vector<vector<double>> I = identity(geo.dimension());
  vector<double> grad_P0(geo.dimension()*geo.elements()[0].size());
  int nbElements = geo.elements()[0].size();
  int nbNodesPerElement = geo.dimension()+1;
  Vec u0_seq = getSequential(u);

  for (int idElement = 0; idElement < geo.elements()[0].size(); idElement++){ /* loop on boundary elements */

    /* Boundary element */
    vector<int> simplex = geo.elements()[0][idElement];
    vector<vector<double>> coordinates(nbNodesPerElement);
    for (int i = 0; i < nbNodesPerElement; i++){
      coordinates[i] = geo.coordinates()[simplex[i]];
    }
    fe.setCoordinates(coordinates);

    /* Get elemental velocity */
    vector<double> u_element(nbNodesPerElement); 
    u_element = geo.getNodalValues(u0_seq, simplex, 1);

    /* Compute grad u = T * u, with T = JT^-1 * grad_phi */
    vector<vector<double>> T = transpose(fe.jacobianInverse()) * G;
    vector<double> grad_u =  T * u_element;

    for (int i = 0; i < geo.dimension(); i++){
      grad_P0[idElement + i*nbElements] = grad_u[i];
    }
  }
  return m_interpolator.interpolateP0_P1(grad_P0);
}

/* Compute P1 gradient on the boundary. For vector fields, return \grad u \cdot n, instead */
Vec Calculus::gradientOnBoundary(Vec u, int bdLabel){
  if (par.nbDofsPerNode()[0] > 1){
    return gradientOnBoundaryForVectorField(u, bdLabel);
  } else {
    return gradientOnBoundaryForScalarField(u, bdLabel);
  }
}

/* Computes \nabla(u) \cdot n */
Vec Calculus::gradientOnBoundaryForVectorField(Vec u, int bdLabel){

  vector<double> gradient_P0(par.nbDofsPerNode()[0]*geo.nbElementsBoundary(), 0.0);
  Vec u0_seq = getSequential(u);

  for (int idElement = 0; idElement < geo.elementsBD(bdLabel).size(); idElement++){ /* loop on boundary elements */
    vector <int> simplexBD = geo.elementsBD(bdLabel)[idElement];
    /* get finite element coordinates */
    vector<vector<double>> coordinatesBD(m_dimension);
    for (int i = 0; i < m_dimension; i++){
      coordinatesBD[i] = geo.coordinates()[simplexBD[i]];
    }
    feBD.setCoordinates(coordinatesBD);

    /* Work with nodal values from tetra associated with current triangle (needed for P1)  */
    vector<int> simplex = geo.getElementFromElementBD(simplexBD); 

    vector<vector<double>> coordinates(m_dimension + 1);
    for (int i = 0; i < m_dimension + 1; i++){
      coordinates[i] = geo.coordinates()[simplex[i]];
    }
    fe.setCoordinates(coordinates);

    vector<double> u_el = geo.getNodalValues(u0_seq, simplex, par.nbDofsPerNode()[0]);

    /* Compute P0 gradient */
    vector<vector<double>> grad_u(m_dimension);
    for (int i = 0; i < m_dimension; i++){
      for (int j = 0; j < par.nbDofsPerNode()[0]; j++){
        grad_u[i].resize(par.nbDofsPerNode()[0], 0.0);
      }
    }
    for (int dof_i = 0; dof_i < m_dimension + 1; dof_i++){
      for (int comp = 0; comp < m_dimension; comp++){
        for (int i = 0; i < par.nbDofsPerNode()[0]; i++){
          grad_u[i][comp] += u_el[par.nbDofsPerNode()[0]*dof_i + i] * dot(getRow(fe.jacobianInverse(), comp), fe.grad_phi(dof_i));
        }
      }
    }
    grad_u = transpose(grad_u);
    feBD.computeNormal();
    vector<double> shearElemental = grad_u*feBD.normal();

    /* Assemble P0 vector */
    int idLabelOffset = 0;
    for (int i = 0; i < findSomething(geo.bdLabels(), bdLabel); i++){
      idLabelOffset += geo.elementsBD()[i].size(); 
    }
    for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
      gradient_P0[par.nbDofsPerNode()[0]*(idLabelOffset + idElement) + comp] = shearElemental[comp];
    }
  }

  /* Interpolate to P1 */
  return m_interpolator.interpolateP0_P1_boundary(gradient_P0);
}

Vec Calculus::gradientOnBoundaryForScalarField(Vec u, int bdLabel){
  if (m_world_rank == 0) cout << "gradient on boundary for scalar field not yet implemented." << endl;
  exit(1);
}

Vec Calculus::divergence(Vec u){
  if (m_world_rank == 0) cout << "CAL: Computing divergence." << endl;
  
  vector<double> divergence_P0(geo.elements()[0].size(), 0.0);

  MasterElement fe; 
  fe.initialize(par, geo.dimension());
  int nbNodesPerElement = geo.dimension()+1;
  Vec u_seq = getSequential(u);
  for (int elementId=0; elementId < geo.elements()[0].size(); elementId++){
    /* Set finite element */
    vector<int> element = geo.elements()[0][elementId];
    vector<vector<double>> coord(nbNodesPerElement);
    for (int i=0; i<nbNodesPerElement; i++){
      coord[i] = geo.coordinates()[element[i]];
    }
    fe.setCoordinates(coord);
    /* Loop on element nodes  */
    for (int dof_i = 0; dof_i < nbNodesPerElement; dof_i++){
      vector<double> u_node(par.nbDofsPerNode()[0]); 
      code = VecGetValues(u_seq, par.nbDofsPerNode()[0], &range(par.nbDofsPerNode()[0]*element[dof_i], par.nbDofsPerNode()[0]*element[dof_i]+par.nbDofsPerNode()[0])[0], &u_node[0]); CHKERR(code);
      divergence_P0[elementId] += dot(u_node, (transpose(fe.jacobianInverse())*fe.grad_phi(dof_i)) );
    }
  }
  VecDestroy(&u_seq);
  return m_interpolator.interpolateP0_P1(divergence_P0);
}

Vec Calculus::wallShearStress(vector<Vec> u, int bdLabel){

  vector<vector<double>> G = transpose(fe.grad_phi());
  vector<vector<double>> I = identity(geo.dimension());
  vector<double> wss(geo.nbElementsBoundary(), 0.0); 
  int nbNodesPerElement = geo.dimension()+1;
  vector<Vec> u0_seq(geo.dimension());
  for (int j = 0; j < geo.dimension(); j++){
    u0_seq[j] = getSequential(u[j]);
  }
//vector<double> asd(geo.elements()[0].size(), 1.0);
//return m_interpolator.interpolateP0_P1(asd);

  vector<vector<double>> normals = transpose(bd.normalsP0comp());
  for (int idElement = 0; idElement < geo.elementsBD(bdLabel).size(); idElement++){ /* loop on boundary elements */

    /* Boundary element */
    vector<int> simplexBD = geo.elementsBD(bdLabel)[idElement];
    vector<vector<double>> coordinatesBD(m_dimension);
    for (int i = 0; i < m_dimension; i++){
      coordinatesBD[i] = geo.coordinates()[simplexBD[i]];
    }
    feBD.setCoordinates(coordinatesBD);
//    feBD.computeNormal();

    /* Inner element */
    vector<int> simplex = geo.getElementFromElementBD(simplexBD); 
    vector<vector<double>> coord(nbNodesPerElement);
    for (int idNode = 0; idNode < nbNodesPerElement; idNode++){
      coord[idNode] = geo.coordinates()[simplex[idNode]];
    }
    fe.setCoordinates(coord);

    /* Compute T = JT^-1 * grad_phi */
    vector<vector<double>> T = transpose(fe.jacobianInverse()) * G;

    /* Get elemental velocity */
    vector<vector<double>> u_element(geo.dimension()); 
    for (int j = 0; j < geo.dimension(); j++){
      u_element[j] = geo.getNodalValues(u0_seq[j], simplex, 1);
    }
    vector<vector<double>> grad_u = T * transpose(u_element);

    int idLabelOffset = 0;
    for (int i = 0; i < findSomething(geo.bdLabels(), bdLabel); i++){
      idLabelOffset += geo.elementsBD()[i].size(); 
    }

    double eta = 0;
    if (par.power_law_n() != 1.0){
      double gamma_dot = sqrt(0.5*trace(grad_u*transpose(grad_u)));
      if (gamma_dot != 0){
        eta = par.power_law_m() * pow(gamma_dot, par.power_law_n() - 1.0);
      } 
      if (eta < par.power_law_m()/10000.0){
        eta = par.power_law_m()/10000.0;
      } else if (eta > 10000.0*par.power_law_m()){
        eta = par.power_law_m()*10000.0;
      }
    } else {
      eta = par.viscosity();
    }

    vector<double> normal = normals[idElement+idLabelOffset];
    for (int i = 0; i < geo.dimension(); i++)
    wss[idElement + idLabelOffset] = eta * 0.5 * norm((I - transpose(transpose(normal))*transpose(normal)) * ( grad_u + transpose(grad_u) ) * normal);
  }
  return m_interpolator.interpolateP0_P1_boundary(wss);
}

vector<Vec> Calculus::split(Vec u){
  vector<Vec> splited(par.nbVariables());
  Vec u0_seq = getSequential(u);
  int m, n;
  code = VecGetOwnershipRange(u, &m, &n); CHKERR(code);
  for (int i = 0; i < par.nbVariables(); i++){
    int offset = 0;
    for (int j = 1; j < i + 1; j++){
      offset += par.nbDofsPerNode()[j-1]*geo.nbVertices;
    }
    int nbDofs = par.nbDofsPerNode()[i]*geo.nbVertices;
    double * sol = new double[nbDofs];
    code = VecGetValues(u0_seq, nbDofs, &range(offset, offset + nbDofs)[0], sol); CHKERR(code);
    vec(splited[i], nbDofs);
    for (int dofId = 0; dofId < nbDofs; dofId++){
      vecSetInsert(splited[i], dofId, sol[dofId]);
    }
    delete [] sol;
    VecAssemblyBegin(splited[i]);
    VecAssemblyEnd(splited[i]);
  }
  VecDestroy(&u0_seq);
  return splited; 
}

vector<Vec> Calculus::decompose_vector_field(Vec u, string mode){
  vector<Vec> u_comps(geo.dimension());
  if (mode == "contiguous"){
    for (int i = 0; i < geo.dimension(); i++){
      u_comps[i] = vec(m_nbVertices);
    }
    Vec u0_seq = getSequential(u);
    for (int i = 0; i < geo.coordinates().size(); i++){
      vector<double> u_node(geo.dimension());
      code = VecGetValues(u0_seq, geo.dimension(), &range(geo.dimension()*i,geo.dimension()*i+geo.dimension())[0], &u_node[0]); CHKERR(code);
      for (int j = 0; j < geo.dimension(); j++){
        vecSetInsert(u_comps[j], i, u_node[j]);
      }
    }
    for (int j = 0; j < geo.dimension(); j++){
      VecAssemblyBegin(u_comps[j]);
      VecAssemblyEnd(u_comps[j]);
    }
  } else {
    if (m_world_rank == 0) cout << "ERROR: Only contiguous mode implemented for function decompose_vector_field." << endl;
    exit(1);
  }
  return u_comps;

}

Vec Calculus::recompose_vector_field(vector<Vec> u, string mode){ /* decompose_vector_field^{-1} */
  int nbDofs = 0;
  assert(u.size() == geo.dimension() && "ERROR: stl vector size in function recompose_vector_field must equal geo.dimension().");
  vector<vector<double>> u_stl(geo.dimension());
  for (int i = 0; i < geo.dimension(); i++){
    int vecSize;
    code = VecGetSize(u[i], &vecSize); CHKERR(code);
    nbDofs = nbDofs + vecSize;
    u_stl[i] = stl(u[i]);
  }
  Vec v = vec(nbDofs);
  if (mode == "contiguous"){
    for (int i = 0; i < geo.nbVertices; i++){
      for (int comp = 0; comp < geo.dimension(); comp++){
        vecSetInsert(v, geo.dimension()*i+comp, u_stl[comp][i]);
      }
    } 
  } else if (mode == "split"){
    if (m_world_rank == 0) cout << "ERROR: split mode not yet implemented for function recompose_vector_field." << endl;
    exit(1);
  } else {
    exit(1);
  }
  VecAssemblyBegin(v);
  VecAssemblyEnd(v);
  return v;
}

Vec Calculus::magnitude(Vec u){

  vector<double> u_stl = stl(u);
  Vec mag = zeros(geo.nbVertices);
  for (int i = 0; i < geo.nbVertices; i++){
    double value = 0.0;
    for (int j = 0; j < par.nbDofsPerNode()[0]; j++){
      value += u_stl[par.nbDofsPerNode()[0]*i + j] * u_stl[par.nbDofsPerNode()[0]*i + j];
    }
    vecSetInsert(mag, i, sqrt(value));
  }
  VecAssemblyBegin(mag);
  VecAssemblyEnd(mag);

  return mag;
}

Vec Calculus::join(vector<Vec> u, IO & inputOutput){

  IO joint_io;
  joint_io = inputOutput;

  int nbDofsTotal = 0;
  vector<int> nbDofVar(par.nbVariables());
  for (int i = 0; i < par.nbVariables(); i++){
    nbDofsTotal += par.nbDofsPerNode()[i]*joint_io.nbVertices();
    nbDofVar[i] = par.nbDofsPerNode()[i]*joint_io.nbVertices();
  }

  Vec joint = vec(nbDofsTotal);
  for (int i = 0; i < par.nbVariables(); i++){
    int offset = 0;
    for (int j = 0; j < i; j++){
      offset += nbDofVar[j];
    }
    Vec u0_seq = getSequential(u[i]);
    int nbDofs = par.nbDofsPerNode()[i]*joint_io.nbVertices();
    double * sol = new double[nbDofs];
    code = VecGetValues(u0_seq, nbDofs, &range(nbDofs)[0], sol); CHKERR(code);
    for (int dofId = 0; dofId < nbDofs; dofId++){
      vecSetInsert(joint, offset + dofId, sol[dofId]);
    }
    delete [] sol;
    VecDestroy(&u0_seq);
  }

  VecAssemblyBegin(joint);
  VecAssemblyEnd(joint);
  return joint; 
}

double Calculus::boundaryIntegralScalar(Vec u, int bdLabel){
  return bd.integral(u, bdLabel);
}
