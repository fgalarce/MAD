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

#include<learning.hpp>

void Learning::initialize(Parameters parameters, Geometry geometry){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "Learning: Initializing learning tools." << endl;
  par = parameters;
  m_geo = geometry;
  m_nbVertices = geometry.nbVertices;
  m_nbDofs = geometry.nbVertices*3;
}

void Learning::initialize(Parameters parameters){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "Learning: Initializing learning tools." << endl;
  par = parameters;
  mkdir(par.dirResults());
}

void Learning::voxelize(){
  if (m_world_rank == 0) cout << "Learning: Voxelizing geometry. " << endl;
  m_voxels = importdata(par.cartesian_grid()); 
  m_nbVoxels = m_voxels.size();
  m_voxelized_geometry.resize(m_nbVoxels);
  this->get_clouds_per_voxel();
  if (par.voxelization_type() == "points"){
    this->voxelize_points();
  } else if (par.voxelization_type() == "volume"){
    this->voxelize_volume();
  }
}

void Learning::get_clouds_per_voxel(){
  m_clouds.resize(m_nbVoxels);
  for (int idVoxel = 0; idVoxel < m_nbVoxels; idVoxel++){
    for (int idPoint = 0; idPoint < m_geo.coordinates().size(); idPoint++){
      vector<double> e1(3, 0.0); e1[0] = 1.0;
      vector<double> e2(3, 0.0); e2[1] = 1.0;
      vector<double> e3(3, 0.0); e3[2] = 1.0;
      if ( m_geo.isInside(m_geo.coordinates()[idPoint], m_voxels[idVoxel], e1, e2, e3)){
        m_clouds[idVoxel].push_back(idPoint);
      }
    }
  }
}

void Learning::voxelize_points(){
  for (int i = 0; i < m_nbVoxels; i++){
    m_voxelized_geometry[i] = m_clouds[i].size(); 
  }
}

void Learning::voxelize_volume(){

//  for (int i = 0; i < m_clouds.size(); i++){
//
//    /* Gather all tetrahedras inside voxel */
//    vector<int> T;
//    for (int j = 0; j < m_clouds[i].size(); j++){ 
//      vector<int> S = m_geo.star(m_clouds[i][j]);
//      for (int k : S){
//        bool push_or_not_push = true;
//        for (int l : T){
//          if (k == l){
//            push_or_not_push = false;
//          }
//        }
//        if (push_or_not_push){
//          T.push_back(k);
//        }
//      }
//    }
//
//    MasterTetrahedron tet;
//    tet.initialize(par); 
//
//    m_voxelized_geometry[i] = 0.0;
//
//    /* Compute and collect volumes */
//    for (int j : T){
//
//      vector<int> tetra = m_geo.tetrahedron()[0][j];
//
//      /* get finite element coordinates */
//      vector<vector<double>> coordinates(4);
//      coordinates[0] = m_geo.coordinates()[tetra[0]];
//      coordinates[1] = m_geo.coordinates()[tetra[1]];
//      coordinates[2] = m_geo.coordinates()[tetra[2]];
//      coordinates[3] = m_geo.coordinates()[tetra[3]];
//
//      
//      tet.setCoordinates(coordinates);
//      m_voxelized_geometry[i] += tet.detJacobian(); 
//    }
//  }
}

void Learning::write_geo_unstructured_grid(){

  string saveName = par.dirResults() + "/g_.vtk";
  if (m_world_rank == 0) cout << "Learning: Writing " << saveName << endl;
  ofstream outfile(saveName.c_str());
  int nx = 6;
  int ny = 50;
  int nz = 6;

  outfile << "# vtk DataFile Version 2.0 " << endl;
  outfile << "Data animation" << endl;
  outfile << "ASCII" << endl;
  outfile << "DATASET STRUCTURED_POINTS" << endl;
  outfile << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
  outfile << "ASPECT_RATIO 1.0 1.0 1.0" << endl;
  outfile << "ORIGIN 0.0 0.0 0.0" << endl;
  outfile << "POINT_DATA " << m_nbVoxels << endl;
  outfile << "SCALARS Field float" << endl;
  outfile << "LOOKUP_TABLE default" << endl;
  for(int h=0; h<m_nbVoxels; h++){
    double value = m_voxelized_geometry[h];
    if(fabs(value) < 1.0e-9){
      value = 0.0;
    }
    outfile << value << endl;
  }
  outfile.close();
}

/* receives squared dissimilarity matrix */
Mat & Learning::MDS(Mat D){
  
  if (m_world_rank == 0) cout << "Learning: Computing multi-dimensional scaling." << endl;

  int nbSamples;
  MatGetSize(D, &nbSamples, NULL);

  if (m_world_rank == 0) cout << "Learning: Assembling X'X = - H D2 H." << endl;
  Mat H = eye(nbSamples, "dense");
  Vec uno = ones(nbSamples);
  code = MatAXPY(H, -1.0/nbSamples, outer(uno), DIFFERENT_NONZERO_PATTERN); CHKERR(code);

  Mat XX = mat(nbSamples, nbSamples, "dense");
  Mat XX0 = mat(nbSamples, nbSamples, "dense");
  code = MatMatMult(H, D, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &XX0); CHKERR(code);
  code = MatMatMult(XX0, H, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &XX); CHKERR(code);
  code = MatScale(XX, -1.0/2.0); CHKERR(code);

  if (m_world_rank == 0) cout << "Learning: Computing low-dimensional representation." << endl;
  EPS eps;
  int nconv;
  EPSCreate(PETSC_COMM_WORLD, &eps);
  EPSSetOperators(eps, XX, NULL);
  EPSSetProblemType(eps, EPS_HEP);
  EPSSetType(eps, EPSKRYLOVSCHUR);
  EPSSetDimensions(eps, nbSamples, PETSC_DECIDE, PETSC_DECIDE);
  EPSSetFromOptions(eps);
  EPSSolve(eps);
  EPSGetConverged(eps, &nconv);
  if (m_world_rank == 0) cout << "EPS: nconv : " << nconv << endl; 

  vector<double> eValues(nbSamples);
  vector<double> sValues;
  int dimension = par.nbModes();
  vector<Vec> V(dimension);

  int k = 0;
  for (int i = 0; i < dimension; i++){

    code = EPSGetEigenvalue(eps, i, &eValues[i], NULL); CHKERR(code);
    if (m_world_rank == 0) cout << "Learning: MDS. V" << i << ". EV = " << eValues[i] << endl;

    /* Compute mode */
    if (eValues[i] > 0){
      vec(V[k], nbSamples);
      code = EPSGetEigenvector(eps, i, V[k], NULL); CHKERR(code);
      sValues.push_back(sqrt(eValues[i]));
      k = k + 1;
    } else {
      V.erase(V.begin() + k);
    }
  }

  Mat U = buildProjector(V);
  Mat Y = mat(sValues.size(), dimension, "dense");
  code = MatMatMult(diag(sValues), transpose(U), MAT_INITIAL_MATRIX, PETSC_DEFAULT, &Y); CHKERR(code);

  saveMat(Y, par.dirResults() + "/Y.m");
//  exportData(par.dirResults() + "/D_hat.m", dissimilarityMatrix(stl(Y)));
  exportData(par.dirResults() + "/eigenValues.txt", eValues);

  return Y;
}

vector<vector<double>> Learning::dissimilarityMatrix(vector<vector<double>> Y){
  int m = Y[0].size();
  vector<vector<double>> D(m);
  for (int i = 0; i < m; i++){
    D[i].resize(m);
    for (int j = 0; j < m; j++){
      D[i][j] = norm(Y[i] - Y[j]);
    }
  }
  return D;
}

/* Place geometry voxelizations on the columns of a matrix */
void Learning::assembleG(){
  if (m_world_rank == 0) cout << "Learning: Assembling matrix from voxelized geometry data. " << endl;

  int m_nbSamples = par.nbMeshes() * (1 + par.noiseIterations());

  int nbVoxels = (importdata1D(par.geometryData() + "/v" + wildcard(0) + "/g.txt")).size(); 
  mat(m_G, nbVoxels, m_nbSamples, "dense");
  int k = 0;
  for (int i = 0; i < par.nbMeshes(); i++){
    vector<double> g = importdata1D(par.geometryData() + "/v" + wildcard(i) + "/g.txt"); 
    MatSetValues(m_G, nbVoxels, &range(nbVoxels)[0], 1, &k, &g[0], INSERT_VALUES);
    if (par.noiseIterations() > 0 ){
      for (int noiseIt = 0; noiseIt < par.noiseIterations(); noiseIt++){
        vector<double> g = importdata1D(par.geometryData() + "/v" + wildcard(i) + "/g" + wildcard(noiseIt) + ".txt"); 
        MatSetValues(m_G, nbVoxels, &range(nbVoxels)[0], 1, &k, &g[0], INSERT_VALUES);
        k++;
      }
    } else{
      k++;
    }
  }
  MatAssemblyBegin(m_G, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m_G, MAT_FINAL_ASSEMBLY);

}

void Learning::assemble_normal_equations(){
  m_pinvG = bestAATinverse(m_G);
//  m_pinvG = bestATAinverse(m_G);
}

vector<double> Learning::compute_weights(int idRow){
  if (m_world_rank == 0) cout << "Learning: Computing weights. " << endl;
  vector<vector<double>> E = importdata(par.errorMatrix());
  MatGetSize(m_G, &m_nbVoxels, NULL);
  vec(m_w, m_nbVoxels);
  code = MatMult(m_pinvG, petsc(E[idRow]), m_w); CHKERR(code);
  Vec res = vec(par.nbMeshes()); 
  Mat tG = mat(par.nbMeshes(), m_nbVoxels);
  code = MatTranspose(m_G, MAT_INITIAL_MATRIX, &tG); CHKERR(code);
  MatMult(tG, m_w, res); 
  VecAXPY(res, -1.0, petsc(E[idRow]));
  VecNorm(res, NORM_2, &m_residual);
  cout << "Learning: Residual = " << m_residual << endl;
  return stl(m_w);
}
