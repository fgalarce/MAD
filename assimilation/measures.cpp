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

#include<measures.hpp>

vector<double> noslip(vector<double> x, double t, Parameters par){
  vector<double> bc(par.nbDofsPerNode()[0], 0.0);
  return bc;
}

void Measures::finalize(){
  if (m_world_rank == 0) cout << "MEMORY: Releasing Measures memory." << endl;
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  for (int i = 0; i < m_rieszRepresenters.size(); i++){
    code = VecDestroy(&m_rieszRepresenters[i]); CHKERRV(code);
  }
}

void Measures::initialize(const Parameters parameters, Geometry & geometry, LinearAlgebra & la, IO & inputOutput){

  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);

  if (m_world_rank == 0) cout << "Measures: Initializing." << endl;
  par = parameters;
  geo = geometry;
  m_la = la;
  m_ip = m_la.ip;
  io = inputOutput;
  /* Compute problem size */
  for (int i = 0; i < par.nbVariables(); i++){
    m_nbDofs = m_nbDofs + geo.nbVertices*par.nbDofsPerNode()[i];
  }

  /* get finite elements */
  if (geo.dimension() == 2){
    elements = geo.triangles();
  } else {
    elements = geo.tetrahedron();
  }

  /* Import measures from file */ 
  vector<vector<double>> image = importdata(par.measures());

  m_imageFOR.resize(3);
  m_imageFOR[0] = par.beamDir();
  m_imageFOR[1] = par.transDir();
  m_imageFOR[2] = cross(m_imageFOR[0], m_imageFOR[1])/norm(cross(m_imageFOR[0], m_imageFOR[1]));

  m_nbVoxels = image.size();
  m_coordinates.resize(m_nbVoxels);

  if (par.modeMeasures() == "CFI" || par.modeMeasures() == "scalarAverage"){
    m_nbComponentsPerVoxel = 1;
  } else if (par.modeMeasures() == "VFI" || par.modeMeasures() == "MRE"){
    m_nbComponentsPerVoxel = 2;
  } else if (par.modeMeasures() == "4d-flow" || par.modeMeasures() == "vectorAverage"){
    m_nbComponentsPerVoxel = 3;
  }
  m_nbMeasures = m_nbComponentsPerVoxel*m_nbVoxels;
  m_measures.resize(m_nbMeasures);

  if (m_world_rank == 0) cout << "Measures: Raw number of voxels: " << m_nbVoxels << endl;
  if (m_world_rank == 0) cout << "Measures: Raw number of measures: " << m_nbMeasures << endl;

  for (int i = 0 ; i < m_nbVoxels; i++){
    m_coordinates[i].resize(geo.dimension());
    for (int comp = 0; comp < geo.dimension(); comp++){
      m_coordinates[i][comp] = image[i][comp];
    }
    for (int comp = 0; comp < m_nbComponentsPerVoxel; comp++){
      m_measures[m_nbComponentsPerVoxel*i + comp] = image[i][geo.dimension() + comp]; 
    }
  }
  m_rieszRepresenters.resize(m_nbMeasures);
  m_rieszRepresentersL2.resize(m_nbMeasures);
  if (par.innerProduct() == "H1" && par.rieszRepresentersFolder() == " "){ 
    /* compute and save L2 representation anyway so the data can be shown on the mesh */ 
    this->computeRieszRepresentersL2();
    m_la.orthonormalize(m_rieszRepresenters); /* use mGS */
    m_rieszRepresentersL2.resize(m_rieszRepresenters.size());
    for (int i = 0; i < m_nbMeasures; i++){
      vec(m_rieszRepresentersL2[i], m_nbDofs);
      VecCopy(m_rieszRepresenters[i], m_rieszRepresentersL2[i]);
      VecDestroy(&m_rieszRepresenters[i]);
    }
    /* re set and re iniliatize representers for H1 space */
    m_rieszRepresenters.resize(m_nbComponentsPerVoxel*m_nbVoxels);
    this->computeRieszRepresentersH1();
    m_la.orthonormalize(m_rieszRepresenters); /* use mGS */
  } else if ((par.innerProduct() == "L2" or par.innerProduct() == "l2") && par.rieszRepresentersFolder() == " ") {
    this->computeRieszRepresentersL2();
    m_la.orthonormalize(m_rieszRepresenters); /* use mGS */
    m_rieszRepresentersL2.resize(m_rieszRepresenters.size());
    for (int i = 0; i < m_nbMeasures; i++){
      vec(m_rieszRepresentersL2[i], m_nbDofs);
      VecCopy(m_rieszRepresenters[i], m_rieszRepresentersL2[i]); /* TODO: avoid variable redunancy in L2 framework */
    }
  } else {
    assert(par.nbMeasures() > 0);
    m_nbMeasures = par.nbMeasures();
    m_rieszRepresenters.resize(m_nbMeasures);
    for (int i = 0; i < m_nbMeasures; i++){
      vec(m_rieszRepresenters[i], m_nbDofs);
      io.loadVector(m_rieszRepresenters[i], par.rieszRepresentersFolder() + "/observer." + wildcard(i) + ".bin");
    }
  }

  if (m_world_rank == 0) cout << "Measures: updated number of measures is " << m_nbMeasures << endl;

  if (!par.useSyntheticMeasures()){
    vec(m_l, m_nbMeasures);
    exit(1);
    code = VecSetValues(m_l, m_nbMeasures, &range(m_nbMeasures)[0], &m_measures[0], INSERT_VALUES); CHKERRV(code);
    code = VecAssemblyBegin(m_l); CHKERRV(code);
    code = VecAssemblyEnd(m_l); CHKERRV(code); 
  }

  if (par.modeMeasures() == "CFI"){
    computeMeasuresMax();
  }
  if  (par.gaussianNoiseLevel() > 0){
    this->computeNoiseIntensity();
  }
}

void Measures::computeMeasuresMax(){
  if (m_world_rank == 0) cout << "Measures: computing maximal measure magnitude for noise model." << endl;
  vector<Vec> rr(m_nbMeasures);
  if (par.nbVariables() > 1){
    int nbDofsFirst = par.nbDofsPerNode()[0]*geo.nbVertices;
    for (int i = 0; i < m_nbMeasures; i++){
      vec(rr[i], nbDofsFirst);
      vector<double> rr_array(nbDofsFirst);
      code = VecGetValues(getSequential(m_rieszRepresenters[i]), nbDofsFirst, &(range(nbDofsFirst))[0], &rr_array[0]); CHKERR(code);
      int rr_low, rr_high;
      code = VecGetOwnershipRange(rr[i], &rr_low, &rr_high); CHKERR(code);
      for(int j = 0; j < nbDofsFirst; j++){
        if (j >= rr_low and j < rr_high){
          VecSetValue(rr[i], j, rr_array[j], INSERT_VALUES);
        }
      }
      code = VecAssemblyBegin(rr[i]); CHKERR(code);
      code = VecAssemblyEnd(rr[i]); CHKERR(code);
    }

  } else {
    for (int i = 0; i < m_nbMeasures; i++){
      vec(rr[i], par.nbDofsPerNode()[0]*geo.nbVertices);
      VecZeroEntries(rr[i]);
      VecAXPY(rr[i], 1.0, m_rieszRepresenters[i]);
    }
  }
  for (int i = par.start(); i < par.end(); i++){
    double li_max = 0.0;
    Vec target = vec(par.nbDofsPerNode()[0]*geo.nbVertices);
//    loadVec(target, par.dirSyntheticField() + "/" + par.variableName()[0] + "." + wildcard(i) + ".vct");
    loadVec(target, par.dirSyntheticField() + "/" + par.variableName()[0] + "." + wildcard(i) + ".scl");
    for (int i = 0; i < m_rieszRepresenters.size(); i++){
      double li = m_ip(target, rr[i], 0); CHKERR(code);
      if (li > li_max){
        li_max = li;
      }
    }
    m_measureMaxValue.push_back(li_max);
  }
}

void Measures::computeNoiseIntensity(){
  m_std_dev = max(m_measureMaxValue) / par.gaussianNoiseLevel();
  if (m_world_rank == 0) cout << "Measures: noise standard deviation = " << m_std_dev << endl;
}

void Measures::computeRieszRepresentersH1(){
  if (m_world_rank == 0) cout << "Measures: Computing H1 Riesz representers li = <omega, u>. \n";
  Boundary boundary;
  boundary.initialize(par, geo);
  
  PetscErrorCode code;

  MasterElement fe; 
  fe.initialize(par, geo.dimension());

  Mat LHS = mat(m_nbDofs, m_nbDofs);
  MatDuplicate(m_ip.matrix, MAT_COPY_VALUES, &LHS);
  MatAssemblyBegin(LHS, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(LHS, MAT_FINAL_ASSEMBLY);

  for (int i = 0; i < geo.bdLabels().size(); i++){
    boundary.Dirichlet(geo.bdLabels()[i], noslip);
  }
  boundary.block(LHS);
  KSP ksp = configureKSP(LHS, "preonly", "lu");
  vector<int> nzRiesz;
  for (int idVoxel = 0; idVoxel < m_nbVoxels; idVoxel++){

    if (m_world_rank == 0) cout << "  Voxel " << wildcard(idVoxel) << ". Center: ";
    for (int compVox = 0; compVox < geo.dimension(); compVox++){
      if (m_world_rank == 0) cout << scientific << m_coordinates[idVoxel][compVox] << " ";
    }

    int coordId = 0;
    vector<int> voxelSupport; voxelSupport.resize(0);
    vector<int> pointsInside;

    for (auto it = geo.coordinates().cbegin(); it != geo.coordinates().cend(); it++){
      if ( geo.isInside(*it, m_coordinates[idVoxel], m_imageFOR[0], m_imageFOR[1], m_imageFOR[2])) {
        for (int comp = 0; comp < geo.dimension(); comp++){
          voxelSupport.push_back(geo.dimension()*coordId+comp);
        }
        pointsInside.push_back(coordId);
      }
      coordId++;
    }

    /* Assemble Riesz representers. Assemble only if the measure will be not discarded afterwards */ 
    if (voxelSupport.size() > 0){
      if (m_world_rank == 0) cout << "Points inside: " << pointsInside.size() << endl;

      /* LHS computation L_i = int_{idVoxel} phi_i \cdot b */
      vector<Vec> L(m_nbComponentsPerVoxel);
      for (int i = 0; i < m_nbComponentsPerVoxel; i++){
        L[i] = zeros(m_nbDofs);
      }
      int low,high;
      code = VecGetOwnershipRange(L[0], &low, &high); CHKERR(code);

      for (int i = 0; i < pointsInside.size(); i++){ /* loop over points inside image voxel */
        for (int comp = 0; comp < geo.dimension(); comp++){
          if (geo.dimension()*i+comp >= low && geo.dimension()*i+comp < high){
            vector<vector<int>> neighbours;
            neighbours = geo.getElementNeighbourhood(pointsInside[i]);
            vector<int> tetraSet = neighbours[0];
            vector<int> indexes = neighbours[1];
            vector<int> partIds = neighbours[2];

            /* Compute int_voxel \phi_i \cdot b */
            double Li;
            for (int j = 0; j < tetraSet.size(); j++){ /* Loop over tetrahedras */
              vector<int> tetra = elements[partIds[j]][tetraSet[j]];
              vector<vector<double>> coordinates(geo.dimension()+1);
              for (int vertex = 0; vertex < geo.dimension() + 1; vertex++){
                coordinates[vertex] = geo.coordinates()[tetra[vertex]];
              }
              fe.setCoordinates(coordinates);
              Li = 0;
              for (int qp = 0; qp < par.nbQuadraturePoints(); qp++){
                Li = Li + fe.phi(qp, indexes[j]) * fe.weights()[qp] * fe.detJacobian();
              }
            }            
            for (int dim = 0; dim < m_nbComponentsPerVoxel; dim++){
              code = VecSetValue(L[dim], (PetscInt)voxelSupport[geo.dimension()*i+comp], m_imageFOR[dim][comp]*Li, INSERT_VALUES); CHKERRV(code);
            }
          }
        }
      }
      for (int j = 0; j < m_nbComponentsPerVoxel; j++){
        if (m_world_rank == 0) cout << "    Solving <omega_i, v> = int_{Omega_i} v cdot e" << j << endl;
        VecAssemblyBegin(L[j]);
        VecAssemblyEnd(L[j]);
        boundary.block(L[j]);
        vec(m_rieszRepresenters[m_nbComponentsPerVoxel*idVoxel + j], m_nbDofs);
        code = KSPSolve(ksp, L[j], m_rieszRepresenters[m_nbComponentsPerVoxel*idVoxel + j]); CHKERR(code); 
        double rrNorm = sqrt(m_ip(m_rieszRepresenters[m_nbComponentsPerVoxel*idVoxel + j], m_rieszRepresenters[m_nbComponentsPerVoxel*idVoxel + j]));
        if (m_world_rank == 0) cout << "    Norm rr(observer): " << rrNorm << endl;
        /* discard RR if norm is zero */
        if (rrNorm > 0) nzRiesz.push_back(m_nbComponentsPerVoxel*idVoxel + j);
        code = VecDestroy(&L[j]); CHKERR(code);
      }
      pointsInside.clear();
    } else {
      if (m_world_rank == 0) cout << ". Support is empty." << endl;
    }
  }
  /* Update amount of measures */
  m_nbMeasures = nzRiesz.size();
  vector<Vec> rieszRepresentersAux(m_nbMeasures);
  for (int i = 0; i < m_nbMeasures; i++){
    rieszRepresentersAux[i] = zeros(m_nbDofs);
    VecAXPY(rieszRepresentersAux[i], 1.0, m_rieszRepresenters[nzRiesz[i]]);
  }
  m_rieszRepresenters.clear();
  m_rieszRepresenters.resize(m_nbMeasures);
  for (int i = 0; i < m_nbMeasures; i++){
    m_rieszRepresenters[i] = zeros(m_nbDofs);
    VecAXPY(m_rieszRepresenters[i], 1.0, rieszRepresentersAux[i]);
  }
  if (m_world_rank == 0) cout << "Measures: Updated number of measures: " << m_nbMeasures << endl;
  code = KSPDestroy(&ksp); CHKERRV(code); 
  code = MatDestroy(&LHS); CHKERRV(code); 
}

void Measures::computeRieszRepresentersL2(){
  if (m_world_rank == 0) cout << "Measures: computing measures Riesz representers li = <omega, u>. \n";

  Boundary boundary;
  boundary.initialize(par, geo);

  int countRR = 0;
  vector<double> CFI_measures;
  for (int idRiesz = 0; idRiesz < (int)m_nbVoxels; idRiesz++){
    if (m_world_rank == 0) cout << "  Voxel " << wildcard(countRR) <<". Center: "; 
    for (int compId = 0; compId < geo.dimension(); compId++){
      if (m_world_rank == 0) cout << scientific << m_coordinates[idRiesz][compId] << " ";
    }
    
    int coordId = 0;
    vector<int> pointsInside;

    for (auto it = geo.coordinates().cbegin(); it != geo.coordinates().cend(); it++){
      if (geo.isInside(*it, m_coordinates[idRiesz], m_imageFOR[0], m_imageFOR[1], m_imageFOR[2])) {
        pointsInside.push_back(coordId);
      }
      coordId++;
    }

    /* Assemble Riesz representers. Assemble only if the measure will be not discarded afterwards */ 
    if (pointsInside.size() > 0){
      if (m_world_rank == 0) cout << ". Support is not empty." << endl;
      vec(m_rieszRepresenters[countRR], m_nbDofs);
      code = VecZeroEntries(m_rieszRepresenters[countRR]); CHKERR(code);

      int m, n;
      code = VecGetOwnershipRange(m_rieszRepresenters[countRR], &m, &n); CHKERR(code);
      for (int i : pointsInside){
        for (int compId = 0; compId < par.nbDofsPerNode()[0]; compId++){
          if (par.nbDofsPerNode()[0]*i+compId >= m && par.nbDofsPerNode()[0]*i+compId < n){
            code = VecSetValue(m_rieszRepresenters[countRR], par.nbDofsPerNode()[0]*i+compId, m_imageFOR[0][compId], INSERT_VALUES); CHKERRV(code);
          }
        }
      } 
      code = VecAssemblyBegin(m_rieszRepresenters[countRR]); CHKERRV(code);
      code = VecAssemblyEnd(m_rieszRepresenters[countRR]); CHKERRV(code);

      /* Get rrNorms. We need to store them for posterior re-scaling of measures */
      double rrNorm = sqrt(m_ip(m_rieszRepresenters[countRR], m_rieszRepresenters[countRR]));
      m_rrNorms.push_back(rrNorm);
      code = VecScale(m_rieszRepresenters[countRR], 1.0/rrNorm);

      if (m_world_rank == 0) cout << " Norm : " << rrNorm << endl;
      countRR++;
      for (int iComp = 0; iComp < m_nbComponentsPerVoxel; iComp++){
        id_nz_rr.push_back(m_nbComponentsPerVoxel*idRiesz + iComp); 
      }

      if (par.modeMeasures() == "VFI" or par.modeMeasures() == "4d-flow"){
        vec(m_rieszRepresenters[countRR], m_nbDofs);
        code = VecZeroEntries(m_rieszRepresenters[countRR]); CHKERR(code);
        code = VecGetOwnershipRange(m_rieszRepresenters[countRR], &m, &n); CHKERR(code);
        for (int i : pointsInside){
          for (int compId = 0; compId < par.nbDofsPerNode()[0]; compId++){
            if (par.nbDofsPerNode()[0]*i+compId >= m && par.nbDofsPerNode()[0]*i+compId < n){
              code = VecSetValue(m_rieszRepresenters[countRR], par.nbDofsPerNode()[0]*i+compId, m_imageFOR[1][compId], INSERT_VALUES); CHKERRV(code);
            }
          }
        } 
        code = VecAssemblyBegin(m_rieszRepresenters[countRR]); CHKERRV(code);
        code = VecAssemblyEnd(m_rieszRepresenters[countRR]); CHKERRV(code);

        /* Get rrNorms. We need to store them for posterior re-scaling of measures */
        double rrNorm = sqrt(m_ip(m_rieszRepresenters[countRR], m_rieszRepresenters[countRR]));
        m_rrNorms.push_back(rrNorm);
        code = VecScale(m_rieszRepresenters[countRR], 1.0/rrNorm);

        if (m_world_rank == 0) cout << " Norm : " << rrNorm << endl;
        countRR++;
      }

      if (par.modeMeasures() == "4d-flow"){
        vec(m_rieszRepresenters[countRR], m_nbDofs);
        code = VecZeroEntries(m_rieszRepresenters[countRR]); CHKERR(code);
        code = VecGetOwnershipRange(m_rieszRepresenters[countRR], &m, &n); CHKERR(code);
        for (int i : pointsInside){
          for (int compId = 0; compId < par.nbDofsPerNode()[0]; compId++){
            if (par.nbDofsPerNode()[0]*i+compId >= m && par.nbDofsPerNode()[0]*i+compId < n){
              code = VecSetValue(m_rieszRepresenters[countRR], par.nbDofsPerNode()[0]*i+compId, m_imageFOR[2][compId], INSERT_VALUES); CHKERRV(code);
            }
          }
        } 
        code = VecAssemblyBegin(m_rieszRepresenters[countRR]); CHKERRV(code);
        code = VecAssemblyEnd(m_rieszRepresenters[countRR]); CHKERRV(code);

        /* Get rrNorms. We need to store them for posterior re-scaling of measures */
//        double rrNorm = sqrt(m_ip(m_rieszRepresenters[countRR], m_rieszRepresenters[countRR]));
        double rrNorm;
        VecDot(m_rieszRepresenters[countRR], m_rieszRepresenters[countRR], &rrNorm);
        m_rrNorms.push_back(rrNorm);
        code = VecScale(m_rieszRepresenters[countRR], 1.0/rrNorm);

        if (m_world_rank == 0) cout << " Norm : " << rrNorm << endl;
        countRR++;
      }

      if (!par.useSyntheticMeasures()){
        CFI_measures.push_back(m_CFI[idRiesz]);
      }
      pointsInside.clear();

    } else {
      if (m_world_rank == 0) cout << ". Support is empty." << endl;
    }
  }
  /* Update amount of measures */
  m_rieszRepresenters.erase(m_rieszRepresenters.cbegin() + countRR, m_rieszRepresenters.cend());
  m_nbMeasures = m_rieszRepresenters.size();
  m_nbVoxels = m_nbMeasures / m_nbComponentsPerVoxel;
}

void Measures::computeSynthetic(int iteration){
  if (m_world_rank == 0) cout << "Measures: computing synthetic measures. " << endl;

  /* Import synthetic field from file */
  vector<double> li(m_rieszRepresenters.size());
  Vec u = vec(m_nbDofs);

  vector<string> to_load(par.nbVariables());
  for (int var = 0; var < par.nbVariables(); var++){
    string end_str;
    if (par.nbDofsPerNode()[var] == 1){
      to_load[var] = par.dirSyntheticField() + "/" + par.variableName()[var] + "." + wildcard(iteration) + ".scl"; 
    } else {
      to_load[var] = par.dirSyntheticField() + "/" + par.variableName()[var] + "." + wildcard(iteration) + ".vct"; 
    }
  }
  io.loadVector(u, to_load);
//  if (par.nbVariables() == 1){
//    io.loadVector(u, par.dirSyntheticField() + "/" + par.variableName()[0] + "." + wildcard(iteration) + ".vct");
//  } else if (par.nbVariables() == 2){
//    io.loadVector(u,  par.dirSyntheticField() + "/" + par.variableName()[0] + "." + wildcard(iteration) + ".vct",
//                      par.dirSyntheticField() + "/" + par.variableName()[1] + "." + wildcard(iteration) + ".scl");
//  } else {
//    cout << "WRONG";
//    exit(1);
//  }

  /* Compute W^T u  */
  for (int i = 0; i < m_rieszRepresenters.size(); i++){
    li[i] = m_ip(u, m_rieszRepresenters[i]);
  }

  random_device device_random_;
  default_random_engine generator_(device_random_());
  normal_distribution<> noise(0.0, m_std_dev);
  for (int i = 0; i < m_nbMeasures; i++){
    li[i] = li[i] + noise(generator_);
  }

  vec(m_l, m_rieszRepresenters.size());
  int low,high;
  code = VecGetOwnershipRange(m_l, &low, &high); CHKERR(code);
  for (int i = 0; i < m_nbMeasures; i++){ 
    if (par.modeMeasures() == "CFI"){
      //li[i] = li[i] - (li[i] / max(m_measureMaxValue)) * li[i] ;
    }
    if (i >= low && i < high){
      code = VecSetValue(m_l, i, li[i], INSERT_VALUES); CHKERRV(code);
    }
  }
  PetscErrorCode code;
  code = VecAssemblyBegin(m_l); CHKERRV(code);
  code = VecAssemblyEnd(m_l); CHKERRV(code);
  code = VecDestroy(&u); CHKERRV(code); 

  /* Also compute measures on the mesh, for review purposes */
  vec(m_measuresOnMesh, m_nbDofs);
  VecZeroEntries(m_measuresOnMesh);
  for (int i = 0; i < m_nbMeasures; i++){
    VecAXPY(m_measuresOnMesh, stl(m_l)[i], m_rieszRepresenters[i]);
  }
}

void Measures::computeSynthetic(double time){
  int iteration;
  if ( abs(time - floor(time)) > 0.5){
    iteration = floor(time/par.CFDtimeStep()) + 1;
  } else {
    iteration = floor(time/par.CFDtimeStep());
  }
  this->computeSynthetic(iteration);
}

vector<vector<vector<double>>> Measures::matrixImage(){

  vector<double> trunked_measures = stl(m_l);
  vector<double> extended_measures(m_nbComponentsPerVoxel*m_coordinates.size(), 0.0);
  for (int i = 0; i < id_nz_rr.size(); i++){
    extended_measures[id_nz_rr[i]] = trunked_measures[i];  
  }
  
  vector<vector<double>> extended_measures_x(par.nbVoxelsX());
  vector<vector<double>> extended_measures_y(par.nbVoxelsX());
  vector<vector<double>> extended_measures_z(par.nbVoxelsX());

  for (int i = 0; i < par.nbVoxelsX(); i++){
    extended_measures_x[i].resize(par.nbVoxelsY());
    extended_measures_y[i].resize(par.nbVoxelsY());
    extended_measures_z[i].resize(par.nbVoxelsY());
    for (int j = 0; j < par.nbVoxelsY(); j++){
      extended_measures_x[i][j] = extended_measures[m_nbComponentsPerVoxel*i*par.nbVoxelsY()+ m_nbComponentsPerVoxel * j];
      extended_measures_y[i][j] = extended_measures[m_nbComponentsPerVoxel*i*par.nbVoxelsY()+ m_nbComponentsPerVoxel * j+1];
      extended_measures_z[i][j] = extended_measures[m_nbComponentsPerVoxel*i*par.nbVoxelsY()+ m_nbComponentsPerVoxel * j+2];
    }
  }
  vector<vector<vector<double>>> image_2d(m_nbComponentsPerVoxel);
  image_2d[0] = extended_measures_x;
  image_2d[1] = extended_measures_y;
  image_2d[2] = extended_measures_z;
  return image_2d;
}

void Measures::addNoise(double (*noise) (double,Parameters), int variableLabel){
  vector<double> measures = stl(m_l);
  for (int i = 0; i < m_nbMeasures; i++){
    measures[i] = measures[i] + noise(measures[i], par);
  }
  PetscErrorCode code;
  code = VecSetValues(m_l, m_nbMeasures, &range(m_nbMeasures)[0], &measures[0], INSERT_VALUES); CHKERRV(code);
  code = VecAssemblyBegin(m_l); CHKERRV(code);
  code = VecAssemblyEnd(m_l); CHKERRV(code);
}
