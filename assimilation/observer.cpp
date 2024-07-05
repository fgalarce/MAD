/*=============================================================================
  This file is part of the code MAD
  Copyright (C) 2017-2024,
    
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

#include<observer.hpp>

void Observer::initialize(Parameters parameters, Geometry & geometry, InnerProduct & innerProduct, IO & inputOutput, Calculus & calc, LinearAlgebra & linearAlgebra, Boundary & boundary){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);

  par = parameters;
  geo = geometry;
  ip = innerProduct;
  io = inputOutput;
  calculus = calc;
  la = linearAlgebra;
  bd = boundary;

  fe.initialize(par, geo.dimension());

  loadData();

  m_nbDofs = 0;
  for (int i = 0; i < par.nbVariables(); i++){
    m_nbDofs = m_nbDofs + geo.nbVertices*par.nbDofsPerNode()[i];
  }

  if (par.rieszRepresentersFolder() != " "){
    m_rrNorms = importdata1D(par.rieszRepresentersFolder() + "/rrNorms.txt");
    m_rieszRepresenters.resize(par.nbMeasures());
    for (int i = 0; i < par.nbMeasures(); i++){
      vec(m_rieszRepresenters[i], m_nbDofs);
      loadVec(m_rieszRepresenters[i], par.rieszRepresentersFolder() + "/rr." + wildcard(i) + ".bin");
    }
  } else {
    buildObserver();
  }
}

void Observer::exportObserver(){
  exportData(m_rrNorms, par.dirResults() + "/rrNorms.txt");
  for (int i=0; i < m_rieszRepresenters.size(); i++){
    saveVec(m_rieszRepresenters[i], par.dirResults() + "/observer." + wildcard(i) + ".bin");
  }
}

void Observer::loadData(){
  if (m_world_rank == 0) cout << "Observer: filtering measures." << endl;

  m_voxelBasis.resize(3);
  m_voxelBasis[0] = par.beamDir();
  m_voxelBasis[1] = par.transDir();
  m_voxelBasis[2] = cross(m_voxelBasis[0], m_voxelBasis[1])/norm(cross(m_voxelBasis[0], m_voxelBasis[1]));

  data = par.scaleUnits()*importdata(par.measures());
  m_nbVoxels = data.size();

  m_voxelCenters.resize(m_nbVoxels);
  for (int i = 0 ; i < m_nbVoxels; i++){
    m_voxelCenters[i].resize(geo.dimension());
    for (int comp = 0; comp < geo.dimension(); comp++){
      m_voxelCenters[i][comp] = data[i][comp];
    }
  }
  int treshold = 0;
  for (int idVoxel = 0; idVoxel < m_nbVoxels; idVoxel++){
    vector<int> pointsInside;
    int nbPoints = 0;
    for (int i = 0; i < geo.coordinates().size(); i++){
      if (geo.isInside(geo.coordinates()[i], m_voxelCenters[idVoxel], m_voxelBasis, par.voxelSize())){
        nbPoints = nbPoints + 1;
        pointsInside.push_back(i);
      }
    }
    if (nbPoints > treshold){
      m_pointsInside.push_back(pointsInside);
    } else {
      zeroVoxels.push_back(idVoxel);
    }
  }

  if (m_world_rank == 0) cout << "Observer: raw number of voxels: " << m_nbVoxels << "." << endl;
   
  erase(data, zeroVoxels);

  m_nbVoxels = data.size(); 
  m_measures.resize(m_nbVoxels*par.measureIt().size());
  m_voxelCenters.resize(m_nbVoxels);
  for (int i = 0 ; i < m_nbVoxels; i++){
    m_voxelCenters[i].resize(geo.dimension());
    m_measures[i].resize(par.measureIt().size());
    for (int comp = 0; comp < geo.dimension(); comp++){
      m_voxelCenters[i][comp] = data[i][comp];
    }
    for (int var = 0; var < par.measureIt().size(); var++){
      m_measures[i][var] = data[i][geo.dimension()+var];
    }
  }
  m_nbMeasures = m_measures.size();
  m_rieszRepresenters.resize(m_nbMeasures);
  m_rieszRepresenters_l2.resize(m_nbMeasures);
  m_rrNorms.resize(m_nbMeasures);
  m_rrNorms_l2.resize(m_nbMeasures);
  if (m_world_rank == 0) cout << "Observer: updated number of voxels: " << m_nbVoxels << "." << endl;
}

void Observer::buildObserver(){
  if (m_world_rank == 0) cout << "Observer: computing Riesz representers li = <omega, u>. \n";
  int count_l2 = 0;
  for (int j = 0; j < par.measureIt().size(); j++){
    if (m_world_rank == 0) cout << "Observer: building l2 representer for variable " << par.variableName()[j] << ". \n";
    int offset = 0;
    for (int k = 0; k < j; k++){
      offset = offset + geo.nbVertices*par.nbDofsPerNode()[k];
    }
    build_l2_Observer(j, offset);
    if (par.innerProduct(j) != "l2"){
      if (m_world_rank == 0) cout << "Observer: building " << par.innerProduct(j) << " representer for variable " << par.variableName()[j] << ". \n";
      buildObserverFromIP(j, offset);
    } else {
      count_l2 = count_l2 + 1;
      for (int i = 0; i < m_nbVoxels; i++){
        m_rieszRepresenters[i+m_nbVoxels*j] = vec(m_rieszRepresenters_l2[i+m_nbVoxels*j]);
      }
    }
  }
  /* if all fields are in l2, we do not use mGS */
  if (count_l2 == par.measureIt().size()){
    m_rrNorms = m_rrNorms_l2;
  } else {
    m_rrNorms = la.orthonormalize(m_rieszRepresenters); /* use mGS */
  }
//  exportObserver();
}

void Observer::build_l2_Observer(int idVar, int offset){
  for (int idVoxel = 0; idVoxel < m_nbVoxels; idVoxel++){
    int idRiesz = idVoxel+m_nbVoxels*par.measureIt()[idVar];
    vec(m_rieszRepresenters_l2[idRiesz], m_nbDofs);
    for (int i : m_pointsInside[idVoxel]){
      vecSet(m_rieszRepresenters_l2[idRiesz], offset + i, 1.0); CHKERRV(code);
    }
    code = VecAssemblyBegin(m_rieszRepresenters_l2[idRiesz]); CHKERRV(code);
    code = VecAssemblyEnd(m_rieszRepresenters_l2[idRiesz]); CHKERRV(code);
    m_rrNorms_l2[idRiesz] = sqrt(ip(m_rieszRepresenters_l2[idRiesz], m_rieszRepresenters_l2[idRiesz]));
    VecScale(m_rieszRepresenters_l2[idRiesz], 1.0/m_rrNorms_l2[idRiesz]);
  } 
}


void Observer::buildObserverFromIP(int varId, int offset){
  Mat LHS = mat(ip.ip_mat());
  for (int i = 0; i < geo.bdLabels().size(); i++){
    bd.Dirichlet(geo.bdLabels()[i], noslip_bc, varId);
  }
  bd.block(LHS);
  configureKSP(ksp, par);
  KSPSetOperators(ksp, LHS, LHS);

  for (int idVoxel = 0; idVoxel < m_nbVoxels; idVoxel++){
    Vec rhs = zeros(m_nbDofs);
    for (int i : m_pointsInside[idVoxel]){ /* loop over points inside image voxel */
      vector<vector<int>> neighbours = geo.getElementNeighbourhood(i);
      vector<int> elementSet = neighbours[0];
      vector<int> indexes = neighbours[1];
      /* Compute int_voxel \phi_i */
      double Li = 0.0;
      for (int j = 0; j < elementSet.size(); j++){ /* Loop over tetrahedras */
        vector<int> element = geo.elements()[0][elementSet[j]];
        vector<vector<double>> coordinates(geo.dimension()+1);
        for (int vertex = 0; vertex < geo.dimension() + 1; vertex++){
          coordinates[vertex] = geo.coordinates()[element[vertex]];
        }
        fe.setCoordinates(coordinates);
        for (int qp = 0; qp < par.nbQuadraturePoints(); qp++){
          Li = Li + fe.phi(qp, indexes[j]) * fe.weights()[qp] * fe.detJacobian();
        }
      }            
      vecSet(rhs, offset + i, Li);
    }
    VecAssemblyBegin(rhs);
    VecAssemblyEnd(rhs);
    bd.block(rhs);
    vec(m_rieszRepresenters[idVoxel + par.measureIt()[varId]*m_nbVoxels], m_nbDofs);
    code = KSPSolve(ksp, rhs, m_rieszRepresenters[idVoxel + par.measureIt()[varId]*m_nbVoxels]); CHKERR(code); 
    double norm_l2_rr = norm(m_rieszRepresenters[idVoxel + par.measureIt()[varId]*m_nbVoxels]);
    if (norm_l2_rr == 0){
      errorMessage("buildObserverFromIP", "0-norm Riesz representer. I suggest to increase the threshold for inner points in voxels. Or to re-do it all with a finer mesh.");
    }
    code = VecDestroy(&rhs); CHKERR(code);
  }
  code = KSPDestroy(&ksp); CHKERRV(code); 
  code = MatDestroy(&LHS); CHKERRV(code); 
}

Vec Observer::synthetic(int iteration){
  vec(m_omega, m_nbDofs);
  vec(m_theta, m_nbMeasures);
  VecZeroEntries(m_omega);
  Vec u_gt = vec(m_nbDofs);
  vector<string> to_load(par.nbVariables());
  for (int j = 0; j < par.nbVariables(); j++){
    to_load[j] = par.dirSyntheticField() + "/" + par.variableName()[j] + "." + wildcard(iteration) + ".scl";
  } 
  io.loadVector(u_gt, to_load);

  for (int j : par.measureIt()){
    for (int i = 0; i < m_nbVoxels; i++){
      double theta = ip(u_gt, m_rieszRepresenters[i+j*m_nbVoxels]);
      vecSetInsert(m_theta, i + j*m_nbVoxels, theta);
      code = VecAXPY(m_omega, theta, m_rieszRepresenters[i+j*m_nbVoxels]); CHKERR(code);
    }
  }
  VecAssemblyBegin(m_theta);
  VecAssemblyEnd(m_theta);
  return m_omega;
}

void Observer::buildMeasures(int iteration){
  if (par.useSyntheticMeasures()){
    synthetic(iteration);
  } else {
    zeros(m_theta, m_nbMeasures);
    for (int i = 0; i < m_nbVoxels; i++){
      for (int var = 0; var < par.measureIt().size(); var++){
        vecSet(m_theta, i+var*m_nbVoxels, data[i][geo.dimension()+var]*m_rrNorms[i+var*m_nbVoxels]);
      }
    }
    code = VecAssemblyBegin(m_theta); CHKERR(code);
    code = VecAssemblyEnd(m_theta); CHKERR(code);
  }
}

vector<Vec> Observer::observations(int iteration){
  Vec u_gt = vec(m_nbDofs);
  vector<string> to_load(par.nbVariables());
  for (int j = 0; j < par.nbVariables(); j++){
    to_load[j] = par.dirSyntheticField() + "/" + par.variableName()[j] + "." + wildcard(iteration) + ".scl";
  } 
  io.loadVector(u_gt, to_load);
  vector<Vec> observations(par.measureIt().size());
  Vec omega = zeros(io.nbVertices()*par.nbVariables());
  for (int j : par.measureIt()){
    for (int i = 0; i < m_nbVoxels; i++){
      int idRiesz = i + j*m_nbVoxels;
      if (par.useSyntheticMeasures()){
        code = VecAXPY(omega, ip(u_gt, m_rieszRepresenters_l2[idRiesz]), m_rieszRepresenters_l2[idRiesz]); CHKERR(code);
      } else {
        code = VecAXPY(omega, m_measures[i][j]*m_rrNorms_l2[idRiesz], m_rieszRepresenters_l2[idRiesz]);  CHKERR(code);
      }
    }
  }
  for (int var = 0; var < par.measureIt().size(); var++){
    observations[var] = calculus.split(omega)[var];
  }
  return observations;
}

vector<Vec> Observer::observations(int iteration, string dirField, string sufix){

//  Vec field = vec(m_nbDofs);
//  vector<string> to_load(par.nbVariables());
//  for (int j = 0; j < par.nbVariables(); j++){
//    to_load[j] = dirField + "/" + par.variableName()[j] + sufix + "." + wildcard(iteration) + ".scl";
//  } 
//  io.loadVector(field, to_load);
//  vector<Vec> observations(par.measureIt().size());
//  Vec omega = zeros(io.nbVertices()*par.nbVariables());
//  for (int j : par.measureIt()){
//    for (int i = 0; i < m_nbVoxels; i++){
//      if (par.useSyntheticMeasures()){
//        code = VecAXPY(omega, ip(field, m_rieszRepresenters_l2[i+j*m_nbVoxels]), m_rieszRepresenters_l2[i+j*m_nbVoxels]); CHKERR(code);
//      } 
//    }
//  }
//
//  for (int var = 0; var < par.measureIt().size(); var++){
//    observations[var] = calculus.split(omega)[var];
//  }
//  return observations;

}

//void Observer::loadData(){
//  if (m_world_rank == 0) cout << "Observer: filtering measures." << endl;
//
//  m_voxelBasis.resize(3);
//  m_voxelBasis[0] = par.beamDir();
//  m_voxelBasis[1] = par.transDir();
//  m_voxelBasis[2] = cross(m_voxelBasis[0], m_voxelBasis[1])/norm(cross(m_voxelBasis[0], m_voxelBasis[1]));
//
//  data = par.scaleUnits()*importdata(par.measures());
//  m_nbVoxels = data.size();
//
//  m_voxelCenters.resize(m_nbVoxels);
//  for (int i = 0 ; i < m_nbVoxels; i++){
//    m_voxelCenters[i].resize(geo.dimension());
//    for (int comp = 0; comp < geo.dimension(); comp++){
//      m_voxelCenters[i][comp] = data[i][comp];
//    }
//  }
//  int treshold = 0;
//  vector<int> zeroVoxels;
//  for (int idVoxel = 0; idVoxel < m_nbVoxels; idVoxel++){
//    vector<int> pointsInside;
//    int nbPoints = 0;
//    for (int i = 0; i < geo.coordinates().size(); i++){
//      if (geo.isInside(geo.coordinates()[i], m_voxelCenters[idVoxel], m_voxelBasis, par.voxelSize())){
//        nbPoints = nbPoints + 1;
//        pointsInside.push_back(i);
//      }
//    }
//    if (nbPoints > treshold){
//      m_pointsInside.push_back(pointsInside);
//    } else {
//      zeroVoxels.push_back(idVoxel);
//    }
//  }
//
//  if (m_world_rank == 0) cout << "Observer: raw number of voxels: " << m_nbVoxels << "." << endl;
//   
//  erase(data, zeroVoxels);
//
//  m_nbVoxels = data.size(); 
//  m_measures.resize(m_nbVoxels*par.measureIt().size());
//  m_voxelCenters.resize(m_nbVoxels);
//  for (int i = 0 ; i < m_nbVoxels; i++){
//    m_voxelCenters[i].resize(geo.dimension());
//    m_measures[i].resize(par.measureIt().size());
//    for (int comp = 0; comp < geo.dimension(); comp++){
//      m_voxelCenters[i][comp] = data[i][comp];
//    }
//    for (int var = 0; var < par.measureIt().size(); var++){
//      m_measures[i][var] = data[i][geo.dimension()+var];
//    }
//  }
//  m_nbMeasures = m_measures.size();
//  m_rieszRepresenters.resize(m_nbMeasures);
//  m_rieszRepresenters_l2.resize(m_nbMeasures);
//  m_rrNorms.resize(m_nbMeasures);
//  m_rrNorms_l2.resize(m_nbMeasures);
//  if (m_world_rank == 0) cout << "Observer: updated number of voxels: " << m_nbVoxels << "." << endl;
//}

void Observer::printRieszRepresenters(){
  for (int i = 0; i < m_nbVoxels; i++){
    for (int j = 0; j < par.measureIt().size(); j++){
      io.writeState(calculus.split(m_rieszRepresenters[i + m_nbVoxels*j])[j], "rr_" + par.variableName()[j] + wildcard(i));
    }
  }
}
 
Vec Observer::loadExpected(string expec_folder){
  vector<vector<double>> data_exp = importdata(expec_folder)*par.scaleUnits();
  erase(data_exp, zeroVoxels);
  Vec theta_exp = zeros(m_nbMeasures);
  for (int i = 0; i < m_nbVoxels; i++){
    for (int var = 0; var < par.measureIt().size(); var++){
        //vecSet(m_theta, i+var*m_nbVoxels, data[i][geo.dimension()+var]*m_rrNorms[i+var*m_nbVoxels]);
      vecSet(theta_exp, i+var*m_nbVoxels, data_exp[i][geo.dimension()+var]*m_rrNorms[i+var*m_nbVoxels]);
    }
  }
  code = VecAssemblyBegin(theta_exp); CHKERR(code);
  code = VecAssemblyEnd(theta_exp); CHKERR(code);
  return theta_exp;
}
