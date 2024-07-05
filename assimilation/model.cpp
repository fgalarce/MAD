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

#include <model.hpp>

void Model::finalize(){

  if (m_world_rank == 0) cout << "MEMORY: Releasing model memory." << endl;

  PetscErrorCode code;
  
  for (int i = 0; i < par.nbTimeWindows(); i++){
    for (int j = 0; j < m_nbModesMapping[i]; j++){
      code = VecDestroy(&m_basis[i][j]); CHKERRV(code);
    }
  }
  if (par.type() == "affine"){
    for (size_t i = 0; i < m_uBar.size(); i++)
      code = VecDestroy(&m_uBar[i]); CHKERRV(code);
  }
}

void Model::initialize(Parameters parameters, Geometry & geometry, LinearAlgebra & la, IO & inputOutput){

  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "MODEL: Initializing model. " << endl;
  par = parameters; 
  geo = geometry;
  m_nbVertices = geo.nbVertices;
  m_la = la;
  io = inputOutput;

  /* indexes for petsc routines */
  m_indexes = range(4*m_nbVertices);
  m_indexes_press = range(m_nbVertices);
  m_indexes_press_mix = range(3*m_nbVertices, 4*m_nbVertices);
  m_indexes_vel = range(m_nbVertices*3);

  if (par.nbTimeWindows() > 1){
    initialize_PBDW_pice_wise();
  } else {
    initialize_PBDW_linear();
  }

  if (m_world_rank == 0) cout << "MODEL: Checking orthonormality" << endl;
  for (int i=0; i<m_basis.size(); i++){
    if(!m_la.checkOrthonormality(m_basis[i])){
      m_la.orthonormalize(m_basis[i]);
    }
  }
}

/* OLD initiliazer kept for compatiblity */
void Model::initialize(Parameters parameters, size_t nbVertices){

  if (m_world_rank == 0) cout << "MODEL: Initializing model. " << endl;
  par = parameters; m_nbVertices = nbVertices;

  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);

  /* indexes for petsc routines */
  m_indexes = range(4*m_nbVertices);
  m_indexes_press = range(m_nbVertices);
  m_indexes_press_mix = range(3*m_nbVertices, 4*m_nbVertices);
  m_indexes_vel = range(m_nbVertices*3);

  if (par.nbTimeWindows() > 1){
    initialize_PBDW_pice_wise();
  } else {
    initialize_PBDW_linear();
  }

}

void Model::initialize_PBDW_linear(){
  int nbDofs = 0; 
  for (int i = 0; i < par.nbVariables(); i++){
    nbDofs = nbDofs + geo.nbVertices*par.nbDofsPerNode()[i];
  }
  m_nbModesMapping.push_back(par.nbModes());
  m_basis.resize(1);
  m_basis[0].resize(par.nbModes());

  for (int i = 0; i < par.nbModes(); i++){
    vec(m_basis[0][i], nbDofs);
//    if (par.nbVariables() == 1){
//      io.loadVector(m_basis[0][i], par.dirModel() + "/" + par.variableName()[0] + "_mode." + wildcard(i) + ".vct");
//    } else if (par.nbVariables() == 2){
//      io.loadVector(m_basis[0][i], par.dirModel() + "/" + par.variableName()[0] + "_mode." + wildcard(i) + ".vct",
//                                   par.dirModel() + "/" + par.variableName()[1] + "_mode." + wildcard(i) + ".scl");
//    } else {
//      assert(par.nbVariables() < 3);
//    }

    vector<string> to_load(par.nbVariables());
    for (int var = 0; var < par.nbVariables(); var++){
      string end_str;
      if (par.nbDofsPerNode()[var] == 1){
        to_load[var] = par.dirModel() + "/" + par.variableName()[var] + "_mode." + wildcard(i) + ".scl"; 
      } else {
        to_load[var] = par.dirModel() + "/" + par.variableName()[var] + "_mode." + wildcard(i) + ".vct"; 
      }
    }
    io.loadVector(m_basis[0][i], to_load);
  }
}

void Model::initialize_PBDW_pice_wise(){

  assert(par.HRmax() > 60 and par.HRmax() < 180);
  assert(par.HRmin() > 40 and par.HRmin() < 160);

  /* Locate patient model window */
  double deltaHR = (par.HRmax() - par.HRmin())/par.nbHeartRateWindows();
  for (int i = 0; i < par.nbHeartRateWindows(); i++){
    if (par.heartRate() >= par.HRmin() + deltaHR*i && par.heartRate() < par.HRmin() + deltaHR*(i+1)){
      m_coord_i = i;
      break; 
    }
  }

  if (m_world_rank == 0) cout << "Patient HR window is " << m_coord_i << " : " << par.HRmin() + deltaHR*m_coord_i << " " << par.HRmin() + deltaHR*(m_coord_i+1) << endl;
  if (par.adaptive_nbModes()){
    vector<vector<double>> nbModesMapping_pre = importdata(par.dirModel() + "/" + par.method() + "_" + par.type()  + "_"+ par.innerProduct() + "/mapping_n.txt");
    m_nbModesMapping.resize(nbModesMapping_pre[m_coord_i].size()); 
    for (int i = 0; i < nbModesMapping_pre[m_coord_i].size(); i++){
      m_nbModesMapping[i] = int(nbModesMapping_pre[m_coord_i][i]);
      if (m_world_rank == 0) cout << "n^* (" << i << ") = " << m_nbModesMapping[i] << endl;
    }
  } else {
    m_nbModesMapping.resize(par.nbTimeWindows());
    for (int i = 0; i < par.nbTimeWindows(); i++){
      m_nbModesMapping[i] = par.nbModes();
    }
  }

  PetscInt nbDofs;
  
  /* Compute number of degrees of freedom for mixed reduced model */
  if (par.method() == "pod_mixed")
    nbDofs = 4*m_nbVertices;
  else 
    nbDofs = 3*m_nbVertices;

  /* Load model in HR window for whole cardiac cycle. Notice this is not very memory consuming. Number of time windows are not greater than ~15*/
  m_basis.resize(par.nbTimeWindows());
  m_delta.resize(par.nbTimeWindows());
  for (int i = 0; i < par.nbTimeWindows(); i++){
    m_basis[i].resize(m_nbModesMapping[i]);

    /* Load bounds for constrained least squares with noise */
    if(par.restrainedLeastSquares()){
      if (m_world_rank == 0) cout << "MODEL: loading bounds for constrained least-squares problem with noise." << endl;
      vector<double> fullDelta = importdata1D(par.dirModel() + "/" + par.method() + "_" + par.type()  + "_"+ par.innerProduct()
                                              + "/window_time" + to_string(i) + "_HR" + to_string(m_coord_i)+ "/delta.txt");
      for (int delta_i = 0; delta_i < m_nbModesMapping[i]; delta_i++){
        m_delta[i].push_back(fullDelta[delta_i]);
      }
    }

    if (m_world_rank == 0) cout << "\nMODEL: Loading time window " << i << endl;
    for (int j = 0; j < m_nbModesMapping[i]; j++){
      vec(m_basis[i][j], nbDofs);
      string filename = par.dirModel() + "/" + par.method() + "_" + par.type()  + "_"+ par.innerProduct()+ "/window_time" + to_string(i) + "_HR" + to_string(m_coord_i)  +"/mode";


      if (par.method() == "pod_mixed" && par.modelFormat() == ".vct"){
        Vec vel = vec(3*m_nbVertices);
        Vec press = vec(m_nbVertices);
        loadVec(vel, filename + "_velocity." + wildcard(j) + ".vct");
        loadVec(press, filename + "_pressure." + wildcard(j) + ".scl");

        vector<double> vel_arr(3*m_nbVertices);
        vector<double> press_arr(m_nbVertices);

        VecGetValues(vel, m_indexes_vel.size(), &m_indexes_vel[0], &vel_arr[0]);
        VecGetValues(press, m_indexes_press.size(), &m_indexes_press[0], &press_arr[0]);


        VecSetValues(m_basis[i][j], m_indexes_vel.size(), &m_indexes_vel[0], &vel_arr[0], INSERT_VALUES);
        VecSetValues(m_basis[i][j], m_indexes_press_mix.size(), &m_indexes_press_mix[0], &press_arr[0], INSERT_VALUES);
        VecAssemblyBegin(m_basis[i][j]);
        VecAssemblyEnd(m_basis[i][j]);
        
        VecDestroy(&vel);
        VecDestroy(&press);
      } else {
        loadVec(m_basis[i][j], filename + "." + wildcard(j) + par.modelFormat());
      }

//      loadVec(m_basis[i][j], par.dirModel() + "/" + par.method() + "_" + par.type()  + "/window_time" + to_string(i) + "_HR" + to_string(m_coord_i)  +"/mode." + wildcard(j) + "." + par.modelFormat());
    }
  }

  /* Load snapshots average */
  if (par.type() == "affine"){
    m_uBar.resize(par.nbTimeWindows());

    for (int i = 0; i < par.nbTimeWindows(); i++){
      vec(m_uBar[i], nbDofs);
      loadVec(m_uBar[i], par.dirModel() + "/" + par.method() + "_" + par.type() + "_"+ par.innerProduct() + "/window_time" + to_string(i) + "_HR" + to_string(m_coord_i) + "_"+ par.innerProduct() + "/average." + par.modelFormat());
//      loadVec(m_uBar[i], par.dirModel() + "/" + par.method() + "_" + par.type() + "/window_time" + to_string(i) + "_HR" + to_string(m_coord_i) + "_"+ par.innerProduct() + "/average." + par.modelFormat());
    }
  }    
}

const vector<Vec> & Model::basis(const double time) {

  if (par.nbTimeWindows() > 1) {

    /* Cardiac cycle time is [0,1] */
    double time_norm  = time*par.heartRate()/60;
    double dt = 1.0/par.nbTimeWindows(); 
    m_indexCurrentBasis = -1;

    /* Find coordinate of time window */
    for (int i = 0; i < par.nbTimeWindows(); i++){
      if (time_norm >= dt*i && time_norm < dt*(i+1)){
        m_indexCurrentBasis = i;
        break;
      }
    }

    if (m_indexCurrentBasis == m_indexPreviousBasis){
      m_updateBasis = false;
    } else {
      m_updateBasis = true;
    }

    if (m_indexCurrentBasis == -1){
      if (m_world_rank == 0) cout << "ERROR : required time is outisde model data-base. Extrapolation is not allowed." << endl; 
      exit(1);
    }

    m_indexPreviousBasis = m_indexCurrentBasis;
    return m_basis[m_indexCurrentBasis];
  } else {
    m_indexCurrentBasis = 0;
    return m_basis[0];
  }
}

const Vec Model::average(){
  return m_uBar[m_indexCurrentBasis]; 
}
