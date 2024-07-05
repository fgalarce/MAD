/*=============================================================================
  This file is part of the code MAD
  Multy-physics for mechanicAl engineering and Data assimilation.
  Copyright (C) 2017-2023,
    
     Felipe Galarce

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

#include<io.hpp>

void IO::initialize(Parameters parameters, string dirSurname, string jobSurname){

  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);

  if (m_world_rank == 0) cout << "IO: Initializing." << endl;

  par = parameters;
  m_geometryData = par.geometryData();
  readGeo();

  m_outputDir = par.dirResults() + dirSurname;
  m_jobName = par.patientName() + jobSurname;
  m_geometryName = par.patientName() + ".geo";

  mkdir(m_outputDir); 
  if (!par.minimalOutput()){ 
    writeCase();
    writeGeo();
  }

  if (m_dimension == 2){
    m_nbElementsBoundary = m_nbEdges;
  } else if (m_dimension == 3){
    m_nbElementsBoundary = m_nbTriangles;
  }

}

void IO::initialize(string geometryData, string dirResults){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "IO: Initializing for mesh reading." << endl;
  m_geometryData = geometryData;
  readGeo();

  if (m_dimension == 2){
    m_nbElementsBoundary = m_nbEdges;
  } else if (m_dimension == 3){
    m_nbElementsBoundary = m_nbTriangles;
  }

  m_outputDir = dirResults;
  m_jobName = "secondaryMesh";
  m_geometryName = m_jobName + ".geo";
  m_geometryName = "mesh.geo";

  mkdir(m_outputDir); 
  if (!par.minimalOutput()){ 
    writeCase();
    writeGeo();
  }
}


void IO::initialize(string geometryData, string dirResults, Parameters parameters){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "IO: Initializing for mesh reading." << endl;
  m_geometryData = geometryData;
  readGeo();

  if (m_dimension == 2){
    m_nbElementsBoundary = m_nbEdges;
  } else if (m_dimension == 3){
    m_nbElementsBoundary = m_nbTriangles;
  }

  par = parameters;
  m_outputDir = dirResults;
  m_jobName = par.patientName();
  m_geometryName = par.patientName() + ".geo";

  mkdir(m_outputDir); 
  if (!par.minimalOutput()){ 
    writeCase();
    writeGeo();
  }
}

void IO::writeGeo(){

  if (m_world_rank != 0){
    return;
  }
 
  string geoFileName = m_outputDir + "/" + m_geometryName;
  ofstream geoFile(geoFileName);

  if (geoFile.is_open()){

    if (m_world_rank == 0) cout << "IO: Writing: " << geoFileName << endl;
   
    geoFile << "Geometry file" << endl;
    geoFile << "Geometry file" << endl;
    geoFile << "node id assign" << endl;
    geoFile << "element id assign" << endl;
    geoFile << "coordinates" << endl;
    geoFile << wildcard8(m_nbVertices) << endl;
   
    for (vector<double> coord : m_coordinates){
      for (int i = 0; i < m_dimension; i++){
        if (coord[i] >= 0){
          geoFile << " ";
        }
        geoFile << scientific << setprecision(5) << coord[i];
      }
      if (m_dimension == 2){
        geoFile << scientific << setprecision(5) << " " << 0.0;
      }
      geoFile << endl;
    }
    geoFile << endl;

    int partId = 1; /* labels */ // TODO: correct this and put correct labels

    /* loop on segments parts */
    for (int i = 0; i < m_edges.size(); i++){
      if (i == 0 || (i > 0 && !par.cleanLabels()) ){
        geoFile << "part " << partId << endl;
        geoFile << "part" << partId << endl;
        geoFile << "bar2" << endl;
        if (par.cleanLabels()){
          int nbTetra = 0;
          for (int kk = 0; kk < m_edges.size(); kk++){
            nbTetra = nbTetra + m_edges[kk].size();
          }
          geoFile << wildcard8(nbTetra) << endl; 
        } else {
          geoFile << wildcard8(m_edges[i].size()) << endl; 
        }
        partId++;
      }
      /* loop on segments */  
      for (vector<int> segment : m_edges[i]){
        geoFile << wildcard8(segment[0]) << wildcard8(segment[1]) <<  endl;
      }
      if (!par.cleanLabels()){
        geoFile << endl;
      }
    }


    /* loop on triangle parts */
    if ((m_dimension == 2 && !par.exportSurfaceOnly()) || m_dimension == 3 ){
      for (int i = 0; i < m_triangles.size(); i++){
        if (i == 0 || (i > 0 && !par.cleanLabels()) ){
          geoFile << "part " << partId << endl;
          geoFile << "part" << partId << endl;
          geoFile << "tria3" << endl;
          if (par.cleanLabels()){
            int nbTetra = 0;
            for (int kk = 0; kk < m_triangles.size(); kk++){
              nbTetra = nbTetra + m_triangles[kk].size();
            }
            geoFile << wildcard8(nbTetra) << endl; 
          } else {
            geoFile << wildcard8(m_triangles[i].size()) << endl; 
          }
          partId++;
        }
        /* loop on triangles */  
        for (vector<int> triangle : m_triangles[i]){
          geoFile << wildcard8(triangle[0]) << wildcard8(triangle[1]) << wildcard8(triangle[2]) << endl;
        }
        if (!par.cleanLabels()){
          geoFile << endl;
        }
      }
    }

    if (m_dimension == 3 && !par.exportSurfaceOnly()){
      /* loop on tetrahedron parts */
      for (int i = 0; i < m_tetrahedron.size(); i++){
        if (i == 0 || (i > 0 && !par.cleanLabels())){
          geoFile << "part " << partId << endl;
          geoFile << "part" << partId << endl;
          geoFile << "tetra4" << endl;
          if (par.cleanLabels()){
            int nbTetra = 0;
            for (int kk = 0; kk < m_tetrahedron.size(); kk++){
              nbTetra = nbTetra + m_tetrahedron[kk].size();
            }
            geoFile << wildcard8(nbTetra) << endl; 
          } else {
            geoFile << wildcard8(m_tetrahedron[i].size()) << endl; 
          }
          partId++;
        }
        /* loop on triangles */  
        for (vector<int> tetrahedron: m_tetrahedron[i]){
          geoFile << wildcard8(tetrahedron[0]) << wildcard8(tetrahedron[1]) << wildcard8(tetrahedron[2]) << wildcard8(tetrahedron[3]) << endl;
        }
      }
    }
  } else {
    if (m_world_rank == 0) cout << "ERROR: impossible to open " + geoFileName << endl;
    exit(1);
  }

}

void IO::writeCase(){
  if (m_world_rank != 0){
    return;
  }

  string caseFileName = m_outputDir + "/" + m_jobName + ".case";
  ofstream caseFile(caseFileName);
  
  if (caseFile.is_open()){ 

    if (m_world_rank == 0) cout << "IO: Writing: " << caseFileName << endl;

    caseFile << "FORMAT" << endl;
    caseFile << "type: ensight" << endl;
    caseFile << "GEOMETRY" << endl;
    caseFile << "model: 1 " <<  m_geometryName << endl;
    caseFile.close();

  } else {
    if (m_world_rank == 0) cout << "ERROR: impossible to open " + caseFileName << endl;
    exit(1);
  }
}

void IO::updateCase(string variable, string typeVariable, string perNodeOrPerElement, double time){

  /* Check if variable already belongs to IO object */
  bool variableExist = false;
  for (int i = 0; i < m_variables.size(); i++){
    if (variable == m_variables[i])
      variableExist = true;
  }
  if (!variableExist){
    m_variables.push_back(variable);
    m_variablesType.push_back(typeVariable);
    m_variablesTypeNodeElement.push_back(perNodeOrPerElement);
  }

  ofstream caseFile(m_outputDir + "/" + m_jobName + ".case");
  
  caseFile << "FORMAT" << endl;
  caseFile << "type: ensight" << endl;
  caseFile << "GEOMETRY" << endl;
  caseFile << "model: 1 " <<  m_geometryName << endl;
  caseFile << "VARIABLE" << endl;

  /* Check if the problem is transient */
  if (time != -1.0){

    for (int i = 0; i < m_variables.size(); i++){
      if (m_variablesType[i] == "scl")
        caseFile << "scalar " + m_variablesTypeNodeElement[i] + ": 1 " + m_variables[i] + " " + m_variables[i] + ".*****." + m_variablesType[i] << endl;
      else
        caseFile << "vector " + m_variablesTypeNodeElement[i] + ": 1 " + m_variables[i] + " " + m_variables[i] + ".*****." + m_variablesType[i] << endl;
    }
    /* Check if time already belongs to IO object */
    if (m_time.size() > 0){
      if (time != m_time.back())
        m_time.push_back(time);
    } else {
      m_time.push_back(time);
    }

    caseFile << "TIME" << endl;
    caseFile << "time set: 1" << endl;
    caseFile << "number of steps: " << m_time.size() << endl; 
    caseFile << "filename start number: 0" << endl; 
    caseFile << "filename increment: 1" << endl;
    caseFile << "time values:" << endl;
    for (double t : m_time)
      caseFile << t << endl;
  } else {
    for (int i = 0; i < m_variables.size(); i++){
      if (m_variablesType[i] == "scl")
        caseFile << "scalar " + m_variablesTypeNodeElement[i] +  ": 1 " + m_variables[i] + " " + m_variables[i] + "." + m_variablesType[i] << endl;
      else
        caseFile << "vector " + m_variablesTypeNodeElement[i] +  ": 1 " + m_variables[i] + " " + m_variables[i] + "." + m_variablesType[i] << endl;
    }
  }
  caseFile.close();
}

void IO::writeState(vector<double> u, string variableName, double time){
  writeState(&u[0], variableName, u.size(), time);
}

void IO::writeState(Vec u, double time, const vector<string>& varName){

  /* use default value if any */
  vector<string> vectorName(par.nbVariables());
  if (varName.size() == 0){
    vectorName = par.variableName();
  } else {
    vectorName = varName;
  }

  /* Get vector size */
  int u_size;
  code = VecGetSize(u, &u_size); CHKERRV(code);
  /* scatter to sequential vector */
  Vec u_seq = getSequential(u);
  for (int i = 0; i < par.nbVariables(); i++){
    if (m_world_rank == 0) cout << "IO: writing " << vectorName[i] << endl;
//    double * solArray = new double[m_nbVertices*par.nbDofsPerNode()[i]];
    vector<double> solArray(m_nbVertices*par.nbDofsPerNode()[i]);
    int offset = 0;
    for (int j = 0; j < i; j++){
      offset = offset + par.nbDofsPerNode()[j]*m_nbVertices;
    }
    vector<int> index = range(offset, offset + par.nbDofsPerNode()[i]*m_nbVertices);
    code = VecGetValues(u_seq, par.nbDofsPerNode()[i]*m_nbVertices, &index[0], &solArray[0]); CHKERRV(code); 
    writeState(&solArray[0], vectorName[i], par.nbDofsPerNode()[i]*m_nbVertices, time);
//    delete [] solArray;
    solArray.clear();
    index.clear();
  }
  VecDestroy(&u_seq);
}

void IO::writeState(Vec u, string variableName, double time){
  /* scatter to sequential vector */
  Vec u_seq = getSequential(u);
  /* Get vector size */
  int u_size;
  code = VecGetSize(u_seq, &u_size); CHKERRV(code);
  double * solArray = new double[u_size];
  code = VecGetValues(u_seq, u_size, &range(u_size)[0], solArray); CHKERRV(code); 
  writeState(solArray, variableName, u_size, time);
  code = MPI_Barrier(MPI_COMM_WORLD); CHKERR(code);
  VecDestroy(&u_seq);
  delete [] solArray;
}

void IO::writeState(double * u, string variableName, int vecSize, double time){
  /* Print with first proc only */
  if (m_world_rank != 0 || par.minimalOutput()){
    return;
  }
  /* Check if the field is scalar or vectorial*/
  string typeVariable;
  if (m_nbVertices == vecSize){
    typeVariable = "scl";
  } else if (m_dimension*m_nbVertices == vecSize){
    typeVariable = "vct";
  } else {
    if (m_world_rank == 0) cout << "Error: ensight writeState(...), Vec size: " << vecSize << " not valid." << endl;
    exit(1);
  }

  /* update case file  */
  updateCase(variableName, typeVariable, "per node", time);

  /* adapt file name */
  string stateFileName;
  if (m_time.size() == 0){ /* Stationary problem*/
    stateFileName = m_outputDir + "/" + variableName + "." + typeVariable;
  } else {
    stateFileName = m_outputDir + "/" + variableName + "." + wildcard(m_time.size()-1) + "." + typeVariable;
  }

  ofstream stateFile(stateFileName);

  if (stateFile.is_open()){
    if (m_world_rank == 0) cout << "IO: Writing: " << stateFileName << endl;
    if (m_nbVertices == vecSize){
      stateFile << "Scalar per node" << endl;
    } else if (m_dimension*m_nbVertices == vecSize){
      stateFile << "Vector per node" << endl;
    }
    for (int i = 0; i < vecSize; i++){
      if (u[i] == -0.0){
        u[i] = +0.0;
      }
      if (u[i] >= 0.0){
        stateFile << " ";
      }
      stateFile << scientific << setprecision(5) << u[i];
      if (m_dimension == 2 and typeVariable == "vct"){
        if ((i+1) % 2 == 0 and ((m_dimension + 1)*m_nbVertices != vecSize)){
          stateFile << scientific << setprecision(5) << " " << 0.0;
        }
        if ((i+1) % 4 == 0){
          stateFile << endl;
        }
      } else if (m_dimension == 3 or typeVariable == "scl"){
        if ((i+1) % 6 == 0){
          stateFile << endl;
        }
      }
    }
    stateFile.close();

  } else {
    if (m_world_rank == 0) cout << "ERROR: impossible to open " + stateFileName << endl;
    exit(1);
  }
  return;  
}

//void IO::writeStateP0(vector<double> u, string variableName, vector<int> partId, double time){
//  double u_array[u.size()];
//  for (int i = 0; i < u.size(); i++){
//    u_array[i] = u[i];
//  }
//  writeState(u_array, variableName, partId, u.size(), time);
//}

void IO::writeStateP0(Vec u, string variableName, double time){

  /* Get vector size */
  int u_size;
  code = VecGetSize(u, &u_size); CHKERRV(code);

  vector<vector<int>> elements;
  if (m_dimension == 2){
    elements = m_triangles[0];
  } else {
    elements = m_tetrahedron[0];
  }

  int nbDofs;
  if (u_size == m_dimension*elements.size()){
    nbDofs = m_dimension * elements.size(); // vector per element
  } else {
    nbDofs = elements.size(); // scalar per element
  }

  /* scatter to sequential vector */
  Vec u_seq = getSequential(u);

  double * solArray = new double[nbDofs];
  int offset = 0;
  vector<int> index = range(nbDofs);
  code = VecGetValues(u_seq, nbDofs, &index[0], solArray); CHKERRV(code); 
  writeStateP0(solArray, variableName, nbDofs, time);
}

void IO::writeStateP0(vector<double> u, string variableName, double time){
  writeStateP0(&u[0], variableName, u.size(), time);
}

void IO::writeStateP0(double * u, string variableName, int vecSize, double time){

  if (m_world_rank != 0){
    return;
  }
  cout << "IO: writing " << variableName << endl;

  /* Get element type and decide if the field is scalar or vectorial */
  string elementType, typeVariable;

  if (m_dimension == 3){
    elementType = "tetra4";
    if (vecSize == m_nbTetrahedron){
      typeVariable = "scl";
    } else {
      typeVariable = "vct";
    }
  } else {
    elementType = "tria3";
    if (vecSize == m_nbTriangles){
      typeVariable = "scl";
    } else {
      typeVariable = "vct";
    }
  }

  /* adapt file name */
  string stateFileName;
  if (m_time.size() == 0) /* Stationary problem*/
    stateFileName = m_outputDir + "/" + variableName + "." + typeVariable;
  else 
    stateFileName = m_outputDir + "/" + variableName + "." + wildcard(m_time.size()-1) + "." + typeVariable;

  /* update case file  */
  updateCase(variableName, typeVariable, "per element", time);

  ofstream stateFile(stateFileName);
  if (!stateFile.is_open()){
    if (m_world_rank == 0) cout << "ERROR: impossible to open " + stateFileName << endl;
    exit(1);}

  if (m_world_rank == 0) cout << "IO: Writing: " << stateFileName << endl;

  int nbComp;
  if (typeVariable == "scl"){
    stateFile << "Scalar per element" << endl;
    nbComp = 1;
  } else {
    stateFile << "Vector per element" << endl;
    if (m_dimension == 2){
      nbComp = 2;
    } else {
      nbComp = 3;
    }
  }

  if (elementType == "tria3"){
    stateFile << "part" << wildcard8(m_triangleLabels[0]) << endl << "tria3" << endl; 
    int counter = 0;
    for (int j = 0; j < nbComp*m_triangles[0].size(); j++){
      int u_index;
      u_index = j;
      if (u[u_index] >= 0){
        stateFile << " ";
      }
      stateFile << scientific << setprecision(5) << u[u_index];
      if ((j+1) % 2 == 0){
        stateFile << " " << 0.0;
      }
      if (((counter+1) % 4) == 0){
        stateFile << endl;
      }
      counter++;
    } 
    stateFile << endl << endl;
  }

  if (elementType == "tetra4"){
    stateFile << "part " << wildcard8(m_tetrahedronLabels[0]) << endl << "tetra4" << endl;; 
    for (int j = 0; j < nbComp*m_tetrahedron[0].size(); j++){

      int u_index;
      u_index = j;

      if (u[u_index] >= 0.0){
        stateFile << " ";
      }
      stateFile << scientific << setprecision(5) << u[u_index];
      if ((j+1) % 6 == 0){
        stateFile << endl;
      }
    } 
    stateFile << endl << endl;
  }

  if (elementType == "tria3" && m_dimension == 3){ /* fill tetrahedras with zeros for correct visualization with Paraview */
    for (int id = 0; id < m_tetrahedron.size(); id++){
      stateFile << "part" << wildcard8(m_tetrahedronLabels[id]) << endl << "tetra4" << endl;; 
      for (int j = 0; j < nbComp*m_triangles[id].size(); j++){
        stateFile << " " << scientific << setprecision(5) << 0.0;
        if ((j+1) % 6 == 0){
          stateFile << endl;}
      }
    }

    for (int id = 0; id < m_edges.size(); id++){
      stateFile << "part" << wildcard8(m_edgeLabels[id]) << endl << "bar2" << endl;; 
      for (int j = 0; j < nbComp*m_edges[id].size(); j++){
        stateFile << " " << scientific << setprecision(5) << 0.0;
        if ((j+1) % 6 == 0){
          stateFile << endl;}
      }
    }
  }
  stateFile.close();
}


//void IO::writeStateP0(double * u, string variableName, vector<int> partId, int vecSize, double time){
//
//  /* Get element type and decide if the field is scalar or vectorial */
//  string elementType = getElementType(partId[0]);
//  string typeVariable;
//
//  if (elementType == "tria3"){
//    if (vecSize == m_nbTriangles)
//      typeVariable = "scl";
//    else 
//      typeVariable = "vct";
//  } else if (elementType == "tetra4"){
//    if (vecSize == m_nbTetrahedron)
//      typeVariable = "scl";
//    else
//      typeVariable = "vct";
//  }
//
//  /* adapt file name */
//  string stateFileName;
//  if (m_time.size() == 0) /* Stationary problem*/
//    stateFileName = m_outputDir + "/" + variableName + "." + typeVariable;
//  else 
//    stateFileName = m_outputDir + "/" + variableName + "." + wildcard(m_time.size()-1) + "." + typeVariable;
//
//  /* update case file  */
//  updateCase(variableName, typeVariable, "per element", time);
//
//  ofstream stateFile(stateFileName);
//  if (!stateFile.is_open()){
//    if (m_world_rank == 0) cout << "ERROR: impossible to open " + stateFileName << endl;
//    exit(1);}
//
//  if (m_world_rank == 0) cout << "IO: Writing: " << stateFileName << endl;
//
//  int nbComp;
//  if (typeVariable == "scl"){
//    stateFile << "Scalar per element" << endl;
//    nbComp = 1;
//  } else {
//    stateFile << "Vector per element" << endl;
//    if (m_dimension == 2){
//      nbComp = 2;
//    } else {
//      nbComp = 3;
//    }
//  }
//
//  for (int id = 0; id < partId.size(); id++){
//
//    if (elementType == "tria3"){
//      stateFile << "part" << wildcard8(partId[id]) << endl << "tria3" << endl; 
//      int counter = 0;
//      for (int j = 0; j < nbComp*m_triangles[id].size(); j++){
//
//        int u_index;
//        if (id == 0){
//          u_index = j;
//        } else {
//          u_index = j + nbComp*m_triangles[id-1].size();
//        }
//        if (u[u_index] >= 0){
//          stateFile << " ";
//        }
//        stateFile << scientific << setprecision(5) << u[u_index];
//        if ((j+1) % 2 == 0 && m_dimension == 2){
//          stateFile << " " << 0.0;
//        }
//        if (m_dimension == 2 && ((counter+1) % 4) == 0){
//          stateFile << endl;
//        }
//        if (m_dimension == 3 && ((counter+1) % 6) == 0){
//          stateFile << endl;
//        }
//        counter++;
//      } 
//      stateFile << endl << endl;
//    }
//
//    if (elementType == "tetra4"){
//      stateFile << "part" << wildcard8(partId[id]) << endl << "tetra4" << endl;; 
//      for (int j = 0; j < nbComp*m_tetrahedron[id].size(); j++){
//
//        int u_index;
//        if (id == 0){
//          u_index = j;
//        } else {
//          u_index = j + nbComp*m_triangles[id-1].size();
//        }
//
//        if (u[u_index] >= 0.0){
//          stateFile << " ";
//        }
//        stateFile << scientific << setprecision(5) << u[u_index];
//        if ((j+1) % 6 == 0){
//          stateFile << endl;
//        }
//      } 
//      stateFile << endl << endl;
//    }
//  }
//
//  if (elementType == "tria3"){ /* fill tetrahedras with zeros for correct visualization with Paraview */
//    for (int id = 0; id < m_tetrahedron.size(); id++){
//      stateFile << "part" << wildcard8(m_tetrahedronLabels[id]) << endl << "tetra4" << endl;; 
//      for (int j = 0; j < nbComp*m_triangles[id].size(); j++){
//        stateFile << " " << scientific << setprecision(5) << 0.0;
//        if ((j+1) % 6 == 0){
//          stateFile << endl;}
//      }
//    }
//
//    for (int id = 0; id < m_edges.size(); id++){
//      stateFile << "part" << wildcard8(m_edgeLabels[id]) << endl << "bar2" << endl;; 
//      for (int j = 0; j < nbComp*m_edges[id].size(); j++){
//        stateFile << " " << scientific << setprecision(5) << 0.0;
//        if ((j+1) % 6 == 0){
//          stateFile << endl;}
//      }
//    }
//  }
//  stateFile.close();
//}

void IO::readGeoEnsight(){
  ifstream geoFile(m_geometryData);
  if (m_world_rank == 0) cout << "IO: Reading: " << m_geometryData << endl;
  string line;
  while(getline(geoFile, line)){
    if (line == "coordinates"){
 
      /* get number of points in the mesh */ 
      getline(geoFile, line);
      m_nbVertices = stoi(line);
      if (m_world_rank == 0) cout << "IO: Number of vertices: " << m_nbVertices << endl;
      m_coordinates.resize(m_nbVertices);

      /* get point coordinates */
      string coordinates;
      for (int i = 0; i < m_nbVertices; i++){
        getline(geoFile, coordinates);
        m_coordinates[i].resize(3);
        m_coordinates[i][0] = stod(coordinates.substr( 0, 12));
        m_coordinates[i][1] = stod(coordinates.substr(12, 12));
        m_coordinates[i][2] = stod(coordinates.substr(24, 12));
    }
  } else {

    string partLine = line.substr(0, line.find(" "));

      if (partLine == "part"){
      
        string partId = line.substr(line.find(" ")+1, line.length());  
    
        /* go foward in the file until element type identifier*/ 
        getline(geoFile, line);
        getline(geoFile, line);

        if (line == "bar2"){

          /* Save part label */
          m_edgeLabels.push_back(stoi(partId)); 

          vector<vector<int>> segments;
        
          /* get number of points in the mesh */ 
          getline(geoFile, line);
          m_nbEdges = stoi(line);
          if (m_world_rank == 0) cout << "IO: Number of segments in part " << m_edgeLabels.back() << ": " << m_nbEdges << endl;
          segments.resize(m_nbEdges);

          /* get point coordinates */
          string segment;
          for (int i = 0; i < m_nbEdges; i++){
            getline(geoFile, segment);
            segments[i].resize(2);
            segments[i][0] = stoi(segment.substr(0, 8));
            segments[i][1] = stoi(segment.substr(8, 8));
        }
        m_edges.push_back(segments);

      } else if (line == "tria3"){

          /* Save part label */
          m_triangleLabels.push_back(stoi(partId)); 

          vector<vector<int>> triangles;
        
          /* get number of points in the mesh */ 
          getline(geoFile, line);
          m_nbTriangles = stoi(line);
          if (m_world_rank == 0) cout << "IO: Number of triangles in part " << m_triangleLabels.back() << ": " << m_nbTriangles << endl;
          triangles.resize(m_nbTriangles);

          /* get point coordinates */
          string triangle;
          for (int i = 0; i < m_nbTriangles; i++){

            getline(geoFile, triangle);
            triangles[i].resize(3);
            triangles[i][0] = stoi(triangle.substr(0, 8));
            triangles[i][1] = stoi(triangle.substr(8, 8));
            triangles[i][2] = stoi(triangle.substr(16, 8));

        }
        m_triangles.push_back(triangles);

      } else if (line == "tetra4"){

          /* Save part label */
          m_tetrahedronLabels.push_back(stoi(partId)); 

          vector<vector<int>> tetrahedrons;
        
          /* get number of points in the mesh */ 
          getline(geoFile, line);
          m_nbTetrahedron = stoi(line);
          if (m_world_rank == 0) cout << "IO: Number of tetrahedron in part " << m_tetrahedronLabels.back() << ": " << m_nbTetrahedron << endl;
          tetrahedrons.resize(m_nbTetrahedron);

          /* get point coordinates */
          string tetrahedron;
          for (int i = 0; i < m_nbTetrahedron; i++){
            getline(geoFile, tetrahedron);
            tetrahedrons[i].resize(4);

            tetrahedrons[i][0] = stoi(tetrahedron.substr(0, 8));
            tetrahedrons[i][1] = stoi(tetrahedron.substr(8, 8));
            tetrahedrons[i][2] = stoi(tetrahedron.substr(16, 8));
            tetrahedrons[i][3] = stoi(tetrahedron.substr(24, 8));

          }
          m_tetrahedron.push_back(tetrahedrons);
        }
      }
    }
  }
}

void IO::readGeoEnsightOld(){
  ifstream geoFile(m_geometryData);
  if (m_world_rank == 0) cout << "IO: Reading: " << m_geometryData << endl;
  string line;
  while(getline(geoFile, line)){
    if (line == "coordinates"){
 
      /* get number of points in the mesh */ 
      getline(geoFile, line);
      m_nbVertices = stoi(line);
      if (m_world_rank == 0) cout << "IO: Number of vertices: " << m_nbVertices << endl;
      m_coordinates.resize(m_nbVertices);

      /* get point coordinates */
      string coordinates;
      for (int i = 0; i < m_nbVertices; i++){
        getline(geoFile, coordinates);
        m_coordinates[i].resize(3);

        m_coordinates[i][0] = stod(coordinates.substr(0, coordinates.find(" ")));
        coordinates = coordinates.substr(coordinates.find(" ") + 1, coordinates.length());

        m_coordinates[i][1] = stod(coordinates.substr(0, coordinates.find(" ")));
        coordinates = coordinates.substr(coordinates.find(" ") + 1, coordinates.length());

        m_coordinates[i][2] = stod(coordinates.substr(0, coordinates.find(" ")));
    }

  } else {

//        getline(geoFile, line);

    string partLine = line.substr(0, line.find(" "));

      if (partLine == "part"){
      
        string partId = line.substr(line.find(" ") + 1, line.length());  
    
        /* go foward in the file until element type identifier*/ 
        getline(geoFile, line);
        getline(geoFile, line);

        if (line == "bar2"){

          /* Save part label */
          m_edgeLabels.push_back(stoi(partId)); 

          vector<vector<int>> segments;
        
          /* get number of points in the mesh */ 
          getline(geoFile, line);
          m_nbEdges = stoi(line);
          if (m_world_rank == 0) cout << "IO: Number of segments in part " << m_edgeLabels.back() << ": " << m_nbEdges << endl;
          segments.resize(m_nbEdges);

          /* get point coordinates */
          string segment;
          for (int i = 0; i < m_nbEdges; i++){
            getline(geoFile, segment);
            segments[i].resize(2);
            segments[i][0] = stoi(segment.substr(0, segment.find(" ")));

            segment = segment.substr(segment.find(" ") + 1, segment.length());
            segments[i][1] = stoi(segment.substr(0, segment.find(" ")));

        }
        m_edges.push_back(segments);

      } else if (line == "tria3"){

          /* Save part label */
          m_triangleLabels.push_back(stoi(partId)); 

          vector<vector<int>> triangles;
        
          /* get number of points in the mesh */ 
          getline(geoFile, line);
          m_nbTriangles = stoi(line);
          if (m_world_rank == 0) cout << "IO: Number of triangles in part " + partId + ": " << m_nbTriangles << endl;
          triangles.resize(m_nbTriangles);

          /* get point coordinates */
          string triangle;
          for (int i = 0; i < m_nbTriangles; i++){

            getline(geoFile, triangle);
            triangles[i].resize(3);

            triangles[i][0] = stoi(triangle.substr(0, triangle.find(" ")));
            triangle = triangle.substr(triangle.find(" ") + 1, triangle.length());

            triangles[i][1] = stoi(triangle.substr(0, triangle.find(" ")));
            triangle = triangle.substr(triangle.find(" ") + 1, triangle.length());

            triangles[i][2] = stoi(triangle.substr(0, triangle.find(" ")));

        }
        m_triangles.push_back(triangles);

      } else if (line == "tetra4"){

          /* Save part label */
          m_tetrahedronLabels.push_back(stoi(partId)); 

          vector<vector<int>> tetrahedrons;
        
          /* get number of points in the mesh */ 
          getline(geoFile, line);
          m_nbTetrahedron = stoi(line);
          if (m_world_rank == 0) cout << "IO: Number of tetrahedron in part " + partId + ": " << m_nbTetrahedron << endl;
          tetrahedrons.resize(m_nbTetrahedron);

          /* get point coordinates */
          string tetrahedron;
          for (int i = 0; i < m_nbTetrahedron; i++){
            getline(geoFile, tetrahedron);
            tetrahedrons[i].resize(4);

            tetrahedrons[i][0] = stoi(tetrahedron.substr(0, tetrahedron.find(" ")));
            tetrahedron = tetrahedron.substr(tetrahedron.find(" ") + 1, tetrahedron.length());

            tetrahedrons[i][1] = stoi(tetrahedron.substr(0, tetrahedron.find(" ")));
            tetrahedron = tetrahedron.substr(tetrahedron.find(" ") + 1, tetrahedron.length());

            tetrahedrons[i][2] = stoi(tetrahedron.substr(0, tetrahedron.find(" ")));
            tetrahedron = tetrahedron.substr(tetrahedron.find(" ") + 1, tetrahedron.length());

            tetrahedrons[i][3] = stoi(tetrahedron.substr(0, tetrahedron.find(" ")));

          }
          m_tetrahedron.push_back(tetrahedrons);
        }
      }
    }
  }
}

void IO::readGeo(){
  /* load geometry */
  check_existence(m_geometryData);
  if (get_filename_extension(m_geometryData) == "geo"){
    try {
      readGeoEnsight();
    } catch(...) { /* kept for compatibility with "not well written but still usable" Ensight mesh files */
      if (m_world_rank == 0) cout << "Geo reader failed. Trying with legacy version." << endl;
      readGeoEnsightOld();
    }
  } else if (get_filename_extension(m_geometryData) == "mesh"){
    if (m_world_rank == 0) cout << "IO: Using Inria mesh reader." << endl;
    readGeoInria();
  }
}

void IO::readGeoInria(){

  ifstream geoFile(m_geometryData);
  if (m_world_rank == 0) cout << "IO: Reading: " << m_geometryData << endl;
  string line;
  getline(geoFile, line);
  string meshVersion = line.substr(line.find("MeshVersionFormatted ")+21, line.length());

  if (stod(meshVersion) != 2){
    MDMA_error("Only inria MeshVersionFormated 2 supported");
  }

  while(getline(geoFile, line)){
    if (line.find("Dimension") != string::npos){
      m_dimension = stoi(line.substr(line.find("Dimension ") + 10, line.length()));
    }
  }

  if (m_world_rank == 0) cout << "GEO: mesh dimension = " << m_dimension << endl;

  if (m_dimension == 1) {
    readGeoInria1D();
  } else if (m_dimension == 2) {
    readGeoInria2D();
  } else if (m_dimension == 3) {
    readGeoInria3D();
  } else {
    MDMA_error("Mesh dimension not valid: d = " + to_string(m_dimension));
  }
}

void IO::readGeoInria1D(){
  ifstream geoFile(m_geometryData);
  string line;

  while(getline(geoFile, line)){
    if (line == "Vertices"){
      /* get number of points in the mesh */ 
      getline(geoFile, line);

      m_nbVertices = stoi(line);
      if (m_world_rank == 0) cout << "IO: Number of vertices: " << m_nbVertices << endl;
      m_coordinates.resize(m_nbVertices);

      /* get point coordinates */
      string coordinates;
      for (int i = 0; i < m_nbVertices; i++){
        getline(geoFile, coordinates);

        m_coordinates[i].resize(3);
        m_coordinates[i][0] = stod(coordinates.substr(0, coordinates.find(" ")));
        coordinates = coordinates.substr(coordinates.find(" ") + 1, coordinates.length());

        m_coordinates[i][1] = stod(coordinates.substr(0, coordinates.find(" ")));
        coordinates = coordinates.substr(coordinates.find(" ") + 1, coordinates.length());

        m_coordinates[i][2] = stod(coordinates.substr(0, coordinates.find(" ")));
      }
    } else {

      if (line == "Edges"){

        getline(geoFile, line);
        m_nbEdges = stoi(line);

        if (m_world_rank == 0) cout << "IO: Number of edges: " << m_nbEdges << endl;
  
        string edge;
        vector<vector<int>> edges(m_nbEdges); // + medit label
        vector<int> labelsEdge(m_nbEdges);

        for (int i = 0; i < m_nbEdges; i++){

          getline(geoFile, edge);

          /* Identify subdomain */
          string reversed(edge.rbegin(), edge.rend());

          edges[i].resize(2);

          edges[i][0] = stoi(edge.substr(0, edge.find(" ")));
          edge = edge.substr(edge.find(" ") + 1, edge.length());

          edges[i][1] = stoi(edge.substr(0, edge.find(" ")));
          edge = edge.substr(edge.find(" ") + 1, edge.length());

          labelsEdge[i] = stoi(edge.substr(0, edge.length()));

          /* Look if there are already elements in same subdomain */
          bool new_subdomain = true; 
          for (int j = 0; j < m_edgeLabels.size(); j++){
            if (labelsEdge[i] == m_edgeLabels[j]){
              new_subdomain = false;
            }
          }
          if (new_subdomain){
            m_edgeLabels.push_back(labelsEdge[i]); 
          }
        }
        m_edges.resize(m_edgeLabels.size());

        for (int i = 0; i < edges.size(); i++){
          m_edges[findSomething(m_edgeLabels, labelsEdge[i])].push_back(edges[i]);
        }
      }
    }
  }
}

void IO::readGeoInria2D(){
  ifstream geoFile(m_geometryData);
  string line;

  while(getline(geoFile, line)){
    if (line == "Vertices"){
      /* get number of points in the mesh */ 
      getline(geoFile, line);

      m_nbVertices = stoi(line);
      if (m_world_rank == 0) cout << "IO: Number of vertices: " << m_nbVertices << endl;
      m_coordinates.resize(m_nbVertices);

      /* get point coordinates */
      string coordinates;
      for (int i = 0; i < m_nbVertices; i++){
        getline(geoFile, coordinates);

        m_coordinates[i].resize(2);
        m_coordinates[i][0] = stod(coordinates.substr(0, coordinates.find(" "))) * par.scaleMesh();
        coordinates = coordinates.substr(coordinates.find(" ") + 1, coordinates.length());

        m_coordinates[i][1] = stod(coordinates.substr(0, coordinates.find(" "))) * par.scaleMesh();
      }
    } else {

      if (line == "Edges"){

        getline(geoFile, line);
        m_nbEdges = stoi(line);

        if (m_world_rank == 0) cout << "IO: Number of edges: " << m_nbEdges << endl;
  
        string edge;
        vector<vector<int>> edges(m_nbEdges); // + medit label
        vector<int> labelsEdge(m_nbEdges);

        for (int i = 0; i < m_nbEdges; i++){

          getline(geoFile, edge);

          /* Identify subdomain */
          string reversed(edge.rbegin(), edge.rend());

          edges[i].resize(3);

          edges[i][0] = stoi(edge.substr(0, edge.find(" ")));
          edge = edge.substr(edge.find(" ") + 1, edge.length());

          edges[i][1] = stoi(edge.substr(0, edge.find(" ")));
          edge = edge.substr(edge.find(" ") + 1, edge.length());

          labelsEdge[i] = stoi(edge.substr(0, edge.find(" ")));

          /* Look if there are already elements in same subdomain */
          bool new_subdomain = true; 
          for (int j = 0; j < m_edgeLabels.size(); j++){
            if (labelsEdge[i] == m_edgeLabels[j]){
              new_subdomain = false;
            }
          }
          if (new_subdomain){
            m_edgeLabels.push_back(labelsEdge[i]); 
          }
        }
        m_edges.resize(m_edgeLabels.size());

        for (int i = 0; i < edges.size(); i++){
          m_edges[findSomething(m_edgeLabels, labelsEdge[i])].push_back(edges[i]);
        }
      }

      if (line == "Triangles"){

        getline(geoFile, line);
        m_nbTriangles = stoi(line);

        if (m_world_rank == 0) cout << "IO: Number of triangles: " << m_nbTriangles << endl;
  
        string triangle;
        vector<vector<int>> triangles(m_nbTriangles); // + medit label
        vector<int> labelsTria(m_nbTriangles);

        for (int i = 0; i < m_nbTriangles; i++){

          getline(geoFile, triangle);

          /* Identify subdomain */
          string reversed(triangle.rbegin(), triangle.rend());

          triangles[i].resize(3);

          triangles[i][0] = stoi(triangle.substr(0, triangle.find(" ")));
          triangle = triangle.substr(triangle.find(" ") + 1, triangle.length());

          triangles[i][1] = stoi(triangle.substr(0, triangle.find(" ")));
          triangle = triangle.substr(triangle.find(" ") + 1, triangle.length());

          triangles[i][2] = stoi(triangle.substr(0, triangle.find(" ")));
          triangle = triangle.substr(triangle.find(" ") + 1, triangle.length());

          labelsTria[i] = stoi(triangle.substr(0, triangle.find(" ")));

          /* Look if there are already elements in same subdomain */
          bool new_subdomain = true; 
          for (int j = 0; j < m_triangleLabels.size(); j++){
            if (labelsTria[i] == m_triangleLabels[j]){
              new_subdomain = false;
            }
          }
          if (new_subdomain){
            m_triangleLabels.push_back(labelsTria[i]); 
          }
        }
        m_triangles.resize(m_triangleLabels.size());

        for (int i = 0; i < triangles.size(); i++){
          m_triangles[findSomething(m_triangleLabels, labelsTria[i])].push_back(triangles[i]);
        }
      }

      if (line == "Quadrilaterals"){

        getline(geoFile, line);
        m_nbQuad = stoi(line);

        if (m_world_rank == 0) cout << "IO: Number of triangles: " << m_nbQuad << endl;
  
        string quad;
        vector<vector<int>> quads(m_nbQuad); // + medit label
        vector<int> labelsQuad(m_nbQuad);

        for (int i = 0; i < m_nbQuad; i++){

          getline(geoFile, quad);

          /* Identify subdomain */
          string reversed(quad.rbegin(), quad.rend());

          quads[i].resize(4);

          quads[i][0] = stoi(quad.substr(0, quad.find(" ")));
          quad = quad.substr(quad.find(" ") + 1, quad.length());

          quads[i][1] = stoi(quad.substr(0, quad.find(" ")));
          quad = quad.substr(quad.find(" ") + 1, quad.length());

          quads[i][2] = stoi(quad.substr(0, quad.find(" ")));
          quad = quad.substr(quad.find(" ") + 1, quad.length());

          quads[i][3] = stoi(quad.substr(0, quad.find(" ")));
          quad = quad.substr(quad.find(" ") + 1, quad.length());

          labelsQuad[i] = stoi(quad.substr(0, quad.find(" ")));

          /* Look if there are already elements in same subdomain */
          bool new_subdomain = true; 
          for (int j = 0; j < m_quadLabels.size(); j++){
            if (labelsQuad[i] == m_quadLabels[j]){
              new_subdomain = false;
            }
          }
          if (new_subdomain){
            m_quadLabels.push_back(labelsQuad[i]); 
          }
        }
        m_quads.resize(m_quadLabels.size());

        for (int i = 0; i < quads.size(); i++){
          m_quads[findSomething(m_quadLabels, labelsQuad[i])].push_back(quads[i]);
        }
      }
    }
  }
}

void IO::readGeoInria3D(){

  ifstream geoFile(m_geometryData);
  string line;

  while(getline(geoFile, line)){
    if (line == "Vertices"){
      /* get number of points in the mesh */ 
      getline(geoFile, line);

      m_nbVertices = stoi(line);
      if (m_world_rank == 0) cout << "IO: Number of vertices: " << m_nbVertices << endl;
      m_coordinates.resize(m_nbVertices);
      m_coordinatesLabels.resize(m_nbVertices);

      /* get point coordinates */
      string coordinates;
      for (int i = 0; i < m_nbVertices; i++){
        getline(geoFile, coordinates);

        m_coordinates[i].resize(3);

        m_coordinates[i][0] = stod(coordinates.substr(0, coordinates.find(" "))) * par.scaleMesh();
        coordinates = coordinates.substr(coordinates.find(" ") + 1, coordinates.length());

        m_coordinates[i][1] = stod(coordinates.substr(0, coordinates.find(" "))) * par.scaleMesh();
        coordinates = coordinates.substr(coordinates.find(" ") + 1, coordinates.length());

        m_coordinates[i][2] = stod(coordinates.substr(0, coordinates.find(" "))) * par.scaleMesh();
        coordinates = coordinates.substr(coordinates.find(" ") + 1, coordinates.length());

        m_coordinatesLabels[i] = stod(coordinates);
      }
    } else {

      if (line == "Triangles"){

        getline(geoFile, line);
        m_nbTriangles = stoi(line);

        if (m_world_rank == 0) cout << "IO: Number of triangles: " << m_nbTriangles << endl;
  
        string triangle;
        vector<vector<int>> triangles(m_nbTriangles); // + medit label
        vector<int> labelsTria(m_nbTriangles);

        for (int i = 0; i < m_nbTriangles; i++){

          getline(geoFile, triangle);

          /* Identify subdomain */
          string reversed(triangle.rbegin(), triangle.rend());

          triangles[i].resize(3);

          triangles[i][0] = stoi(triangle.substr(0, triangle.find(" ")));
          triangle = triangle.substr(triangle.find(" ") + 1, triangle.length());

          triangles[i][1] = stoi(triangle.substr(0, triangle.find(" ")));
          triangle = triangle.substr(triangle.find(" ") + 1, triangle.length());

          triangles[i][2] = stoi(triangle.substr(0, triangle.find(" ")));
          triangle = triangle.substr(triangle.find(" ") + 1, triangle.length());

          labelsTria[i] = stoi(triangle.substr(0, triangle.find(" ")));

          /* Look if there are already elements in same subdomain */
          bool new_subdomain = true; 
          for (int j = 0; j < m_triangleLabels.size(); j++){
            if (labelsTria[i] == m_triangleLabels[j]){
              new_subdomain = false;
            }
          }
          if (new_subdomain){
            m_triangleLabels.push_back(labelsTria[i]); 
          }
        }
        m_triangles.resize(m_triangleLabels.size());

        for (int i = 0; i < triangles.size(); i++){
          m_triangles[findSomething(m_triangleLabels, labelsTria[i])].push_back(triangles[i]);
        }
      }
      if (line == "Tetrahedra"){

        getline(geoFile, line);
        m_nbTetrahedron = stoi(line);

        if (m_world_rank == 0) cout << "IO: Number of tetrahedron: " << m_nbTetrahedron << endl;

        string tetra;
        vector<vector<int>> tetrahedron(m_nbTetrahedron); // + medit label
        vector<int> labelsTetra(m_nbTetrahedron);

        for (int i = 0; i < m_nbTetrahedron; i++){

          getline(geoFile, tetra);

          tetrahedron[i].resize(4);

          tetrahedron[i][0] = stoi(tetra.substr(0, tetra.find(" ")));
          tetra = tetra.substr(tetra.find(" ") + 1, tetra.length());

          tetrahedron[i][1] = stoi(tetra.substr(0, tetra.find(" ")));
          tetra = tetra.substr(tetra.find(" ") + 1, tetra.length());

          tetrahedron[i][2] = stoi(tetra.substr(0, tetra.find(" ")));
          tetra = tetra.substr(tetra.find(" ") + 1, tetra.length());

          tetrahedron[i][3] = stoi(tetra.substr(0, tetra.find(" ")));

          labelsTetra[i] = stoi(tetra.substr(tetra.find(" ") + 1, tetra.length()));

          /* Look if there are already elements in same subdomain */
          bool new_subdomain = true; 
          for (int j = 0; j < m_tetrahedronLabels.size(); j++){
            if (labelsTetra[i] == m_tetrahedronLabels[j]){
              new_subdomain = false;
            }
          }
          if (new_subdomain){
            m_tetrahedronLabels.push_back(labelsTetra[i]); 
          }
        }
        m_tetrahedron.resize(m_tetrahedronLabels.size());

        for (int i = 0; i < tetrahedron.size(); i++){
          m_tetrahedron[findSomething(m_tetrahedronLabels, labelsTetra[i])].push_back(tetrahedron[i]);
        }
      }
    }
  }
}

string IO::getElementType(int partId){

//  for (int i = 0; i < m_edgeLabels.size(); i++){
//    if (partId == m_edgeLabels[i])
//      return "bar2";
//  }

  for (int i = 0; i < m_triangleLabels.size(); i++){
    if (partId == m_triangleLabels[i])
      return "tria3";
  }

  for (int i = 0; i < m_tetrahedronLabels.size(); i++){
    if (partId == m_tetrahedronLabels[i])
      return "tetra4";
  }

  return "ELEMENT_TYPE_NON_RECOGNIZED";

}

void IO::writeVTKgrid(vector<double> data, int nx, int ny, int nz, string name,
                      double aspect_ratio_x, double aspect_ratio_y, double aspect_ratio_z){

  string saveName = par.dirResults() + "/" + name + ".vtk";

  if (m_world_rank == 0) cout << "IO: Writing " << saveName << endl;

  ofstream outfile(saveName.c_str());

  outfile << "# vtk DataFile Version 2.0 " << endl;
  outfile << "Data animation" << endl;
  outfile << "ASCII" << endl;
  outfile << "DATASET STRUCTURED_POINTS" << endl;
  outfile << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
  outfile << "ASPECT_RATIO " << aspect_ratio_x << " " << aspect_ratio_y << " " << aspect_ratio_z << endl;
  outfile << "ORIGIN 0.0 0.0 0.0" << endl;
  outfile << "POINT_DATA " << data.size() << endl;
  outfile << "SCALARS Field float" << endl;
  outfile << "LOOKUP_TABLE default" << endl;
  for(int h=0; h<data.size(); h++){
    double value = data[h];
    if(fabs(value) < 1.0e-9){
      value = 0.0;
    }
    outfile << value << endl;
  }
  outfile.close();
}

Vec IO::P1bP1_to_P1P1(Vec u){
  Vec u_seq = getSequential(u);
  double * p1_vec = new double[m_nbVertices*par.nbDofsPerNode()[0]];
  double * p1_scl = new double[m_nbVertices*par.nbDofsPerNode()[1]];
  int nbBubbles = m_nbTriangles*par.nbDofsPerNode()[0];
  code = VecGetValues(u_seq, m_nbVertices*par.nbDofsPerNode()[0], 
                              &range(m_nbVertices*par.nbDofsPerNode()[0])[0], p1_vec); CHKERR(code);
  code = VecGetValues(u_seq, m_nbVertices*par.nbDofsPerNode()[1],
                              &range(m_nbVertices*par.nbDofsPerNode()[0], 
                              nbBubbles + m_nbVertices*par.nbDofsPerNode()[0] + m_nbVertices*par.nbDofsPerNode()[1])[0], p1_scl); CHKERR(code);
  Vec P1_P1 = vec(m_nbVertices*(par.nbDofsPerNode()[0] + par.nbDofsPerNode()[1]));
  VecSetValues(P1_P1, m_nbVertices*par.nbDofsPerNode()[0], 
                      &range(m_nbVertices*par.nbDofsPerNode()[0])[0],
                      p1_vec, INSERT_VALUES);
  VecSetValues(P1_P1, m_nbVertices*par.nbDofsPerNode()[1], 
                      &range(m_nbVertices*par.nbDofsPerNode()[0], 
                             m_nbVertices*par.nbDofsPerNode()[0] + m_nbVertices*par.nbDofsPerNode()[1])[0],
                      p1_scl, INSERT_VALUES);
  return P1_P1;
}

void IO::loadVector(Vec u, vector<string> filename){

  int low, high;
  code = VecGetOwnershipRange(u, &low, &high); CHKERR(code);
  vector<Vec> u_split(par.nbVariables());
  vector<vector<double>> u_stl(par.nbVariables());
  vector<int> nbDofsVar(par.nbVariables());
  int nbDofs = 0;
  for (int i = 0; i < par.nbVariables(); i++){
    nbDofsVar[i] = par.nbDofsPerNode()[i]*m_nbVertices;
    nbDofs += nbDofsVar[i];
    vec(u_split[i], nbDofsVar[i]);
    loadVec(u_split[i], filename[i]);
    u_stl[i] = stl(u_split[i]);

    int offset = 0;
    for (int j = 0; j < i; j++){
      offset += nbDofsVar[j];
    }
    for (int j = 0; j < nbDofsVar[i]; j++){
      if (offset + j >= low && offset + j < high){
        code = VecSetValue(u, offset + j, u_stl[i][j], INSERT_VALUES); CHKERR(code);
      }
    }
  }
  code = VecAssemblyBegin(u); CHKERR(code);
  code = VecAssemblyEnd(u); CHKERR(code);
  double normVector = norm(u);
  if (m_world_rank == 0) cout << "IO: norm loaded vector(s) = " << normVector << endl;
}

void IO::loadVector(Vec u, string filename0, string filename1){

  int nbDofs0 = par.nbDofsPerNode()[0]*m_nbVertices;
  int nbDofs1 = par.nbDofsPerNode()[1]*m_nbVertices;
  int nbDofs = nbDofs0 + nbDofs1;

  /* -- This part of the code is dedicated to ensight ordering of 2d vector as if they were 3d */
  int ensight_vct_numbers0 = 0;
  int ensight_vct_numbers1 = 0;

  if (par.nbDofsPerNode()[0] > 1){
    ensight_vct_numbers0 = ensight_vct_numbers0 + m_nbVertices*3;
  } else {
    ensight_vct_numbers0 = ensight_vct_numbers0 + m_nbVertices;
  }
  if (par.nbDofsPerNode()[1] > 1){
    ensight_vct_numbers1 = ensight_vct_numbers1 + m_nbVertices*3;
  } else {
    ensight_vct_numbers1 = ensight_vct_numbers1 + m_nbVertices;
  }
  /* -- */

  Vec u3D = zeros(ensight_vct_numbers0);
  Vec p3D = zeros(ensight_vct_numbers1);
  loadVec(u3D, filename0);
  loadVec(p3D, filename1);
  vector<double> stl_u_3D = stl(u3D);
  vector<double> stl_p_3D = stl(p3D);
  vector<double> stl_u(nbDofs0);
  vector<double> stl_p(nbDofs1);

  if (par.nbDofsPerNode()[0] > 1){
    for (int i = 0; i < m_nbVertices; i++){
      for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){ /* this jumps the third coordinate always saved in ensight */
        stl_u[par.nbDofsPerNode()[0]*i+comp] = stl_u_3D[3*i+comp];
      }
    }
  } else {
    for (int i = 0; i < m_nbVertices; i++){
      stl_u[i] = stl_u_3D[i];
    }
  }

  if (par.nbDofsPerNode()[1] > 1){
    for (int i = 0; i < m_nbVertices; i++){
      for (int comp = 0; comp < par.nbDofsPerNode()[1]; comp++){ /* this jumps the third coordinate always saved in ensight */
        stl_p[par.nbDofsPerNode()[1]*i+comp] = stl_p_3D[3*i+comp];
      }
    }
  } else {
    for (int i = 0; i < m_nbVertices; i++){
      stl_p[i] = stl_p_3D[i];
    }
  }
  int low, high;
  code = VecGetOwnershipRange(u, &low, &high); CHKERR(code);
  for (int i = 0; i < nbDofs; i++){
    if (i >= low && i < high && i < nbDofs0){
      code = VecSetValue(u, i, stl_u[i], INSERT_VALUES); CHKERR(code);
    }
    if (i >= low && i < high && i > nbDofs0){
      code = VecSetValue(u, i, stl_p[i - nbDofs0], INSERT_VALUES); CHKERR(code);
    }
  }
  code = VecAssemblyBegin(u); CHKERR(code);
  code = VecAssemblyEnd(u); CHKERR(code);
  double normVector = norm(u);
  if (m_world_rank == 0) cout << "IO: norm loaded vector(s) = " << normVector << endl;
}

void IO::loadVector(Vec u, string filename){

  if (this->dimension() == 2){

    int nbDofs = par.nbDofsPerNode()[0]*m_nbVertices;
    /* -- This part of the code is dedicated to ensight ordering of 2d vector as if they were 3d */
    int ensight_vct_numbers;
    if (par.nbDofsPerNode()[0] > 1){
      ensight_vct_numbers = m_nbVertices*3;
    } else {
      ensight_vct_numbers = m_nbVertices;
    }
    /* -- */

    Vec u3D = zeros(ensight_vct_numbers);
    loadVec(u3D, filename);
    vector<double> stl_u_3D = stl(u3D);
    vector<double> stl_u(nbDofs);

    for (int i = 0; i < m_nbVertices; i++){
      for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){ /* this jumps the third coordinate always saved in ensight */
        stl_u[par.nbDofsPerNode()[0]*i+comp] = stl_u_3D[3*i+comp];
      }
    }
    int low, high;
    VecGetOwnershipRange(u, &low, &high);
    for (int i = 0; i < nbDofs; i++){
      if (i >= low && i < high){
        code = VecSetValue(u, i, stl_u[i], INSERT_VALUES); CHKERR(code);
      }
    }
    code = VecAssemblyBegin(u); CHKERR(code);
    code = VecAssemblyEnd(u); CHKERR(code);
  } else {

    if(get_filename_extension(filename) == "scl"){
      int vecSize;
      VecGetSize(u, &vecSize); 
      if (m_nbVertices != vecSize){
        if (m_world_rank == 0) {
          cout << "IO: wrong vector input size when loading: " << filename << endl;
          exit(1);
        }
      }
      loadVec(u, filename);
    } else {
      loadVec(u, filename);
    }
  }
}

//void IO::writeStateP0(Vec u, string variableName, vector<int> partId, double time){
//  /* Get vector size */
//  PetscInt u_size;
//  code = VecGetSize(u, &u_size); CHKERRV(code);
//  PetscInt indexes[u_size];
//  PetscScalar u_array[u_size];
//  /* Build indexes array and get values */
//  for (PetscInt i = 0; i < u_size; i++){
//    indexes[i] = i;
//  }
//  code = VecGetValues(u, u_size, indexes, u_array); CHKERRV(code); 
//  this->writeState(u_array, variableName, partId, u_size, time);
//}
//
