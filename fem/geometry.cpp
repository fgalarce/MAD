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

#include <geometry.hpp>

void Geometry::initialize(const IO & io, vector<double> T){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "GEO: Initializing geometry." << endl;
  m_dimension = io.dimension();
  deployGeo(io); 
  assert(T.size() == m_dimension * m_coordinates.size() );
  for (int i=0; i<m_coordinates.size(); i++){
    for (int comp=0;comp<m_dimension;comp++){
      m_coordinates[i][comp] += T[m_dimension*i+comp];
    }
  }

  if (m_dimension == 1){
    m_elements = m_edges;
    m_elementLabels = m_edgeLabels;
  } else if (m_dimension == 2){
    m_elements = m_triangles;
    m_elementsBD = m_edges;
    m_elementLabels = m_triangleLabels;
    m_elementLabelsBD = m_edgeLabels;
  } else if (m_dimension == 3){
    m_elements = m_tetrahedron;
    m_elementsBD = m_triangles;
    m_elementLabels = m_tetrahedronLabels;
    m_elementLabelsBD = m_triangleLabels;
  }
}

void Geometry::initialize(const IO & io){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "GEO: Initializing geometry." << endl;
  m_dimension = io.dimension();

  deployGeo(io); 

  if (m_dimension == 1){
    m_elements = m_edges;
    m_elementLabels = m_edgeLabels;
  } else if (m_dimension == 2){
    m_elements = m_triangles;
    m_elementsBD = m_edges;
    m_elementLabels = m_triangleLabels;
    m_elementLabelsBD = m_edgeLabels;
  } else if (m_dimension == 3){
    m_elements = m_tetrahedron;
    m_elementsBD = m_triangles;
    m_elementLabels = m_tetrahedronLabels;
    m_elementLabelsBD = m_triangleLabels;
  }
}

void Geometry::initialize(Parameters parameters, const IO & io){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "GEO: Initializing geometry." << endl;
  /* parse parameters */
  par = parameters;
  m_dimension = io.dimension();
  deployGeo(io);
  if (m_dimension == 1){
    m_elements = m_edges;
    m_elementLabels = m_edgeLabels;
  } else if (m_dimension == 2){
    m_elements = m_triangles;
    m_elementsBD = m_edges;
    m_elementLabels = m_triangleLabels;
    m_elementLabelsBD = m_edgeLabels;
  } else if (m_dimension == 3){
    m_elements = m_tetrahedron;
    m_elementsBD = m_triangles;
    m_elementLabels = m_tetrahedronLabels;
    m_elementLabelsBD = m_triangleLabels;
  }
  assembleMass();
  computeMeshVolume();
  print_geo();
}

void Geometry::assembleMass(){

  if (m_world_rank == 0) cout << "GEO: Computing Mass matrix. " << endl;
  mat(m_M, nbVertices, nbVertices); 
  int nbNodesPerElement = m_dimension+1;

  fe.initialize(par, m_dimension);
  for (int partId = 0; partId < m_elements.size(); partId++){
    for (int feId = 0; feId < m_elements[partId].size(); feId++){ /* loop on elements */

      vector<int> simplex = m_elements[partId][feId];
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int i = 0; i < nbNodesPerElement; i++){
        coordinates[i] = m_coordinates[simplex[i]];
      }
      fe.setCoordinates(coordinates);

      if (feId % (m_elements[partId].size() / m_verbose) == 0){
        if (m_world_rank == 0) cout << "  Elementary matrix for element: " << feId << "/" << m_elements[partId].size() - 1 << endl;
      }

      for (int i = 0; i < nbNodesPerElement; i++){
        for (int j = 0; j < nbNodesPerElement; j++){
          double Mij = fe.mass(i,j);
          matSet(m_M, simplex[i], simplex[j], Mij);
        }
      }
    }
  }
  MatAssemblyBegin(m_M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m_M, MAT_FINAL_ASSEMBLY);
  double mNorm = norm(m_M);
  if (m_world_rank == 0) cout << "GEO: Mass matrix norm: " << mNorm << endl;

}

void Geometry::deployGeo(const IO & io){
  if (m_dimension == 1){
    deployGeo1D(io);
  } else if (m_dimension == 2){
    deployGeo2D(io);
  } else if (m_dimension == 3){
    deployGeo3D(io);
  }
}

void Geometry::deployGeo1D(const IO & io){
  m_coordinates = io.coordinates();
  m_coordinatesLabels = io.coordinatesLabels();
  m_edges.resize(io.edges().size());

  /* Change start labeling from 1 (ensight) to 0 */
  for (int i = 0; i < m_edges.size(); i++){ /* loop on ensight parts */
    m_edges[i].resize(io.edges()[i].size());
    m_edges[i] = io.edges()[i] - (int)1;
  }
  nbVertices = io.nbVertices();
}

void Geometry::deployGeo2D(const IO & io){
  m_coordinates = io.coordinates();
  m_coordinatesLabels = io.coordinatesLabels();
  m_triangles.resize(io.triangles().size());
  m_edges.resize(io.edges().size());

  m_edgeLabels = io.edgeLabels();
  m_triangleLabels = io.triangleLabels();

  /* Change start labeling from 1 (ensight) to 0 */
  for (int i = 0; i < m_triangles.size(); i++){ /* loop on ensight parts */
    m_triangles[i].resize(io.triangles()[i].size());
    m_triangles[i] = io.triangles()[i] - (int)1;
  }
  for (int i = 0; i < m_edges.size(); i++){ /* loop on ensight parts */
    m_edges[i].resize(io.edges()[i].size());
    m_edges[i] = io.edges()[i] - (int)1;
  }

  /* Get boundary nodes */
  if (m_world_rank == 0) cout << "GEO: Collecting boundary nodes." << endl;
  m_boundaryNodes.resize(m_edges.size());
  m_nbElementsBoundary = 0;
  for (int i = 0; i < m_edges.size(); i++){ /* loop on ensight parts */
    m_nbElementsBoundary += m_edges[i].size();
    for (vector<int> edge : m_edges[i]){
      for (int j = 0; j < 2; j++){
        vector<int>::iterator it = std::find(m_boundaryNodes[i].begin(), m_boundaryNodes[i].end(), edge[j]);
        if (it != m_boundaryNodes[i].end()){
        } else {
          m_boundaryNodes[i].push_back(edge[j]);
        }
      } 
    }
    if (m_world_rank == 0) cout << "GEO: Number of nodes in boundary part " << m_edgeLabels[i] << ": " << m_boundaryNodes[i].size() << endl;
  }
 
  m_nbBoundaryNodes = 0; 
  for (int i = 0; i < m_boundaryNodes.size(); i++){
    m_nbBoundaryNodes += m_boundaryNodes[i].size();
  }

  m_nbTria = 0;
  for (int i = 0; i < m_triangles.size(); i++){
    m_nbTria += m_triangles[i].size();
    if (m_world_rank == 0) cout << "GEO: Number of triangles in domain part " << m_triangleLabels[i] << ": " << m_triangles[i].size() << endl;
  }
  nbVertices = io.nbVertices();
}

void Geometry::deployGeo3D(const IO & io){
  m_coordinates = io.coordinates();
  m_coordinatesLabels = io.coordinatesLabels();
  m_triangles.resize(io.triangles().size());
  m_tetrahedron.resize(io.tetrahedron().size());
  m_edges.resize(io.edges().size());

  m_edgeLabels = io.edgeLabels();
  m_triangleLabels = io.triangleLabels();
  m_tetrahedronLabels = io.tetrahedronLabels();

  /* Change start labeling from 1 (ensight) to 0 */
  for (int i = 0; i < m_triangles.size(); i++){ /* loop on ensight parts */
    m_triangles[i].resize(io.triangles()[i].size());
    m_triangles[i] = io.triangles()[i] - (int)1;
  }
  for (int i = 0; i < m_edges.size(); i++){ /* loop on ensight parts */
    m_edges[i].resize(io.edges()[i].size());
    m_edges[i] = io.edges()[i] - (int)1;
  }
  for (int i = 0; i < m_tetrahedron.size(); i++){ /* loop on ensight parts */
    m_tetrahedron[i].resize(io.tetrahedron()[i].size());
    for (int j = 0; j < m_tetrahedron[i].size(); j++){
      m_tetrahedron[i][j].resize(io.tetrahedron()[i][j].size());
      for (int k = 0; k < m_tetrahedron[i][j].size(); k++){
        m_tetrahedron[i][j][k] = io.tetrahedron()[i][j][k] - 1;
      }
    }
  }

  /* Get boundary nodes */
  if (m_world_rank == 0) cout << "GEO: Collecting boundary nodes." << endl;
  m_boundaryNodes.resize(m_triangles.size());
  m_nbElementsBoundary = 0;
  for (int i = 0; i < m_triangles.size(); i++){ /* loop on ensight parts */
    m_nbElementsBoundary += m_triangles[i].size();
    for (vector<int> tria : m_triangles[i]){
      for (int j = 0; j < 3; j++){
        vector<int>::iterator it = std::find(m_boundaryNodes[i].begin(), m_boundaryNodes[i].end(), tria[j]);
        if (it != m_boundaryNodes[i].end()){
        } else {
          m_boundaryNodes[i].push_back(tria[j]);
        }
      } 
    }
    if (m_world_rank == 0) cout << "GEO: Number of nodes in boundary part " << m_triangleLabels[i] << ": " << m_boundaryNodes[i].size() << endl;
  }
 
  m_nbBoundaryNodes = 0; 
  for (int i = 0; i < m_boundaryNodes.size(); i++){
    m_nbBoundaryNodes += m_boundaryNodes[i].size();
  }

  m_nbTetrahedra = 0;
  for (int i = 0; i < m_tetrahedron.size(); i++){
    m_nbTetrahedra += m_tetrahedron[i].size();
    if (m_world_rank == 0) cout << "GEO: Number of tetrahedra in domain part " << m_tetrahedronLabels[i] << ": " << m_tetrahedron[i].size() << endl;
  }
  nbVertices = io.nbVertices();
}

bool Geometry::isInside(const vector<double> & point, const vector<double> & voxelCenter, 
                        const vector<vector<double>> & voxelBasis, const vector<double> & voxelSize){
  return isInside(point, voxelCenter, voxelBasis[0], voxelBasis[1], voxelBasis[2], voxelSize[0], voxelSize[1], voxelSize[2]);
} 

bool Geometry::isInside(const vector<double> & point, const vector<double> & voxelCenter, 
                        const vector<vector<double>> & voxelBasis, double voxelSizeX, double voxelSizeY, double voxelSizeZ){
  return isInside(point, voxelCenter, voxelBasis[0], voxelBasis[1], voxelBasis[2], voxelSizeX, voxelSizeY, voxelSizeZ);
} 

bool Geometry::isInside(const vector<double> & point, const vector<double> & voxelCenter, const vector<double> & beamDir, 
                        const vector<double> & transDir, const vector<double> & zDir){ 
  return isInside(point, voxelCenter, beamDir, transDir, zDir, par.voxelSize()[0], par.voxelSize()[1], par.voxelSize()[2]);
}

/* return true if point is inside sample volume centered at vocelCenter*/
bool Geometry::isInside(const vector<double> & point, const vector<double> & voxelCenter, 
                        const vector<double> & beamDir, const vector<double> & transDir, const vector<double> & zDir, 
                        double voxelSizeX, double voxelSizeY, double voxelSizeZ){ 

  double epsL = voxelSizeX/1.0e9;
  double epsT = voxelSizeY/1.0e9;
  double epsZ = voxelSizeZ/1.0e9;  

  vector<double> voxelCenter3d(3,0.0);
  vector<double> point3d(3, 0.0);
  for (int comp = 0; comp < m_dimension; comp++){
    voxelCenter3d[comp] = voxelCenter[comp];
    point3d[comp] = point[comp];
  }

  if (norm( proj( point3d - voxelCenter3d, beamDir )) <= voxelSizeX/2.0-epsL &&
      norm( proj( point3d - voxelCenter3d, transDir)) <= voxelSizeY/2.0-epsT &&
      norm( proj( point3d - voxelCenter3d, zDir )) <= voxelSizeZ/2.0-epsZ ){
    return true;
  } else {
    return false; 
  }
}

/* In d=1: return segment given a point. In d=2: return triangle given a segment. In d=3: return tetrahedra given a triangle. */
vector<int> Geometry::getElementFromElementBD(vector<int> simplexBD){
  vector<int> partId_simplexId;

  for (int i = 0; i < m_elements.size(); i++){ /* loop on labels */
    for (int j = 0; j < m_elements[i].size(); j++){ /* loop on elements */
      int count = 0;
      /* check if current tetrahedron contains all simplexBD points */
      for (int pointId = 0; pointId < m_dimension; pointId++){
        vector<int>::iterator it = std::find(m_elements[i][j].begin(), m_elements[i][j].end(), simplexBD[pointId]);
        if (it != m_elements[i][j].end()){
          count++;
        }
      }
      if (count == m_dimension){
        /* Colect simplex part and label */
        partId_simplexId.push_back(i);
        partId_simplexId.push_back(j);
      }
    }
  }

  vector<int> simplex(m_dimension+1, -1);
  for (int i = 0; i < m_dimension + 1; i++){
    simplex[i] = m_elements[partId_simplexId[0]][partId_simplexId[1]][i];
  }
  
  if (simplex[0] == -1){
    if (m_world_rank == 0) cout << "Tetrahedra not found from triangle with points: " << simplex << ". There is likely a mesh numbering issue. " << endl;
  }
  return simplex;
}

vector<int> Geometry::getTetraFromTriangle(vector<int> triangle){

  if (m_dimension == 2){
    if (m_world_rank == 0) cout << "GEO: function not implemented in 2D yet.\n";
    exit(1);
  }
  vector<int> partId_tetraId;

  for (int i = 0; i < m_tetrahedron.size(); i++){ /* loop on labels */
    for (int j = 0; j < m_tetrahedron[i].size(); j++){ /* loop on elements */
      int count = 0;
      /* check if current tetrahedron contains all triangle points */
      for (int pointId = 0; pointId < 3; pointId++){
        vector<int>::iterator it = std::find(m_tetrahedron[i][j].begin(), m_tetrahedron[i][j].end(), triangle[pointId]);
        if (it != m_tetrahedron[i][j].end()){
          count++;
        }
      }
      if (count == 3){
        /* Colect tetra part and label */
        partId_tetraId.push_back(i);
        partId_tetraId.push_back(j);
        return partId_tetraId;
      }
    }
  }

  if (m_world_rank == 0) cout << "Tetrahedra not found from triangle with points: " << endl;
  if (m_world_rank == 0) cout << triangle << endl;
  return partId_tetraId;
  
}

vector<vector<int>> Geometry::getElementNeighbourhood(int pointId){

  vector<int> neighbours;
  vector<int> index;
  vector<int> partId;

  for (int i = 0; i < m_elements.size(); i++){ /* loop on labels */
    for (int j = 0; j < m_elements[i].size(); j++){ /* loop on elements */
      /* check if current tetrahedron contains pointId */
      vector<int>::iterator it = std::find(m_elements[i][j].begin(), m_elements[i][j].end(), pointId);
      if (it != m_elements[i][j].end()){
        neighbours.push_back(j);
        /* Get position of point in array. This is useful to select the right shape function when integrating during RR routines */
        index.push_back(std::distance(m_elements[i][j].begin(), it));
        /* save part id */
        partId.push_back(i);
      }
    }
  }

  /* Colect neighbours and indexes */
  vector<vector<int>>  neighbours_index_part;
  neighbours_index_part.push_back(neighbours);
  neighbours_index_part.push_back(index);
  neighbours_index_part.push_back(partId);
  return neighbours_index_part;
}

vector<vector<int>> Geometry::getTetraNeighbourhood(int pointId){

  vector<int> neighbours;
  vector<int> index;
  vector<int> partId;
  for (int i = 0; i < m_tetrahedron.size(); i++){ /* loop on labels */
    for (int j = 0; j < m_tetrahedron[i].size(); j++){ /* loop on elements */
      /* check if current tetrahedron contains pointId */
      vector<int>::iterator it = std::find(m_tetrahedron[i][j].begin(), m_tetrahedron[i][j].end(), pointId);
      if (it != m_tetrahedron[i][j].end()){
        neighbours.push_back(j);
        /* Get position of point in array. This is useful to select the right shape function when integrating during RR routines */
        index.push_back(std::distance(m_tetrahedron[i][j].begin(), it));

        /* save part id */
        partId.push_back(m_tetrahedronLabels[i]);
      }
    }
  }

  /* Colect neighbours and indexes */
  vector<vector<int>>  neighbours_index_part;
  neighbours_index_part.push_back(neighbours);
  neighbours_index_part.push_back(index);
  neighbours_index_part.push_back(partId);

  return neighbours_index_part;
}

vector<vector<int>> Geometry::getTriaNeighbourds(int pointId){
  vector<vector<int>> partId_neighbourdId;
  vector<int> partId;
  vector<int> neighbours;
  for (int i = 0; i < m_triangles.size(); i++){ /* loop on labels */
    for (int j = 0; j < m_triangles[i].size(); j++){ /* loop on elements */
      /* check if current triangle contains pointId */
      vector<int>::iterator it = std::find(m_triangles[i][j].begin(), m_triangles[i][j].end(), pointId);
      if (it != m_triangles[i][j].end()){

        /* save part id */
        partId.push_back(m_triangleLabels[i]);
        /* collect neighbour id */
        neighbours.push_back(j);
      }
    }
  }
  partId_neighbourdId.push_back(partId);
  partId_neighbourdId.push_back(neighbours);
  return partId_neighbourdId;
}

vector<vector<int>> Geometry::getTriaNeighbourdsDEPRECATED(int pointId, int i){
  vector<vector<int>> partId_neighbourdId;
  vector<int> partId;
  vector<int> neighbours;
  vector<int> partIDdeprecated;
//  for (int i = 0; i < m_triangles.size(); i++){ /* loop on labels */
    for (int j = 0; j < m_triangles[i].size(); j++){ /* loop on elements  */
      /* check if current triangle contains pointId */
      vector<int>::iterator it = std::find(m_triangles[i][j].begin(), m_triangles[i][j].end(), pointId);
      if (it != m_triangles[i][j].end()){

        /* save part id */
        partId.push_back(m_triangleLabels[i]);
        /* collect neighbour id */
        neighbours.push_back(j);
        /* save part id */
        partIDdeprecated.push_back(i);
      }
    }
//  }
  partId_neighbourdId.push_back(partId);
  partId_neighbourdId.push_back(neighbours);
  partId_neighbourdId.push_back(partIDdeprecated);
  return partId_neighbourdId;
}

vector<int> Geometry::getPointStarBD(int pointId, int i){
  vector<int> neighbours;
  vector<vector<vector<int>>> elements;
  if (m_dimension == 3){
    elements = m_triangles;
  } else if (m_dimension == 2){
    elements = m_edges;
  }
  for (int j = 0; j < elements[i].size(); j++){ /* loop on elements  */
    /* check if current triangle contains pointId */
    vector<int>::iterator it = std::find(elements[i][j].begin(), elements[i][j].end(), pointId);
    if (it != elements[i][j].end()){
      /* collect neighbour id */
      neighbours.push_back(j);
    }
  }
  return neighbours;
}

vector<vector<int>> Geometry::getPointStarBD(int pointId){
  vector<int> neighbours;
  vector<int> partIds;
  vector<vector<vector<int>>> elements;
  for (int i = 0; i < m_elementsBD.size(); i++){
    for (int j = 0; j < m_elementsBD[i].size(); j++){ /* loop on elements  */
      /* check if current triangle contains pointId */
      vector<int>::iterator it = std::find(m_elementsBD[i][j].begin(), m_elementsBD[i][j].end(), pointId);
      if (it != m_elementsBD[i][j].end()){
        /* collect neighbour id */
        neighbours.push_back(j);
        partIds.push_back(i);
      }
    }
  }
  vector<vector<int>> output;
  output.push_back(neighbours);
  output.push_back(partIds);
  return output;
}


void Geometry::computeMeshVolume(){

  Vec uno = zeros(nbVertices);
  code = VecSet(uno, 1.0); CHKERR(code);
  code = VecAssemblyBegin(uno); CHKERR(code);
  code = VecAssemblyEnd(uno);  CHKERR(code);
  Vec M_uno = vec(nbVertices);
  code = MatMult(m_M, uno, M_uno); CHKERR(code);
  code = VecDot(M_uno, uno, &m_meshVolume); CHKERR(code);
}

double Geometry::volume(){

  if (m_world_rank == 0) cout << "GEO: Computing mesh volume.\n";
  m_meshVolume = 0.0;
  for (int partId = 0; partId < m_tetrahedron.size(); partId++) { /* loop on parts */

    if (m_world_rank == 0) cout << "GEO: Part " << partId << endl;
 
    for (vector<int> tetra : m_tetrahedron[partId]){ /* loop on tetra */

      /* get finite element coordinates */
      vector<double> a = m_coordinates[tetra[0]];
      vector<double> b = m_coordinates[tetra[1]];
      vector<double> c = m_coordinates[tetra[2]];
      vector<double> d = m_coordinates[tetra[3]];

      m_meshVolume += 1.0/6.0 * abs(dot(a-d, cross((b-d), (c-d))));
    }
  }

  m_mesh_volume_is_calculated = true;
  return m_meshVolume;
}


double Geometry::computeBoundaryArea(){

  double surfaceArea = 0.0;
  for (int i=0; i<m_triangles.size(); i++){
    surfaceArea += computeBoundaryArea(m_triangleLabels[i]);
  }
  return surfaceArea;
}

double Geometry::computeBoundaryArea(int labelBD){

  if (m_world_rank == 0) cout << "GEO: Computing mesh surface measure part " << labelBD << ".\n";
  double surfaceArea = 0.0;

  for (vector<int> tria : m_triangles[findSomething(m_triangleLabels, labelBD)]){ /* loop on tria */
    vector<double> a = m_coordinates[tria[0]] - m_coordinates[tria[2]];
    vector<double> b = m_coordinates[tria[0]] - m_coordinates[tria[1]];
    surfaceArea += 1.0/2.0 * norm(cross(a, b));
  }
  if (m_world_rank == 0) cout << "GEO: Mesh surface: " << surfaceArea << " cm^2" <<endl;
  return surfaceArea;
}

int Geometry::nearest_neighbourd(vector<double> point){

  /* nearest-neighbourd search */ 
  int NN;
  double distanceST0 = 999999999999; 
  for (int j=0; j < m_coordinates.size(); j++){
    double distanceST1 = norm(m_coordinates[j] - point);
    if (distanceST1 < distanceST0){
      distanceST0 = distanceST1;
      NN = j;
    }
  }

  return NN;
}

vector<int> Geometry::star(int pointId){

  vector<int> star;
  for (int i=0; i<m_tetrahedron[0].size(); i++){ /* loop on elements */
    /* check if current tetrahedron contains pointId */
    for (int j=0; j<m_tetrahedron[0][i].size(); j++){
      if (m_tetrahedron[0][i][j] == pointId){
        star.push_back(i);
      }
    }
  }
  assert(star.size() > 0);
  return star;
}

vector<double> Geometry::getCenter(int labelBD){
  vector<double> center(m_dimension);
  for (int i = 0; i < this->boundaryNodes(labelBD).size(); i++){
    center = center + this->coordinates()[this->boundaryNodes(labelBD)[i]];
  }
  center = center / ((double)this->boundaryNodes(labelBD).size());
  return center;
}

double Geometry::getRadius(int labelBD, vector<double> center){
  double radius = 0.0;
  for (int i = 0; i < this->boundaryNodes(labelBD).size(); i++){
    double candidate_radius = norm(this->coordinates()[this->boundaryNodes(labelBD)[i]] - center);
    if (candidate_radius > radius){
      radius = candidate_radius;
    }
  }
  return radius;
}

vector<double> Geometry::getNodalValues(Vec u0_sequential, vector<int> simplex, int nbDofsPerNode){

  int nbNodesPerElement = simplex.size();

  /* loc2glob dof mapping */
  vector<int> loc2glob(nbNodesPerElement*nbDofsPerNode);
  for (int i = 0; i < nbNodesPerElement; i++){
    for (int comp = 0; comp < nbDofsPerNode; comp++){
      loc2glob[nbDofsPerNode*i+comp] = nbDofsPerNode*simplex[i]+comp;
    }
  }

  vector<double> u_el(nbNodesPerElement*nbDofsPerNode);
  PetscErrorCode code;
  code = VecGetValues(u0_sequential, nbNodesPerElement*nbDofsPerNode, &loc2glob[0], &u_el[0]); CHKERR(code);
  return u_el;
}


void Geometry::print_geo(){
  if (m_dimension == 3){
    for (int i = 0; i < bdLabels().size(); i++){
      double area = computeBoundaryArea(bdLabels()[i]);
      if (m_world_rank == 0) cout << "GEO: boundary area at label " << bdLabels()[i] << ": " << area << endl;  
    }
    computeMeshVolume();
    if (m_world_rank == 0) cout << "GEO: mesh volume: "<< m_meshVolume << endl;  
  }
}

Vec Geometry::boundaryLabelFunction(){
  Vec bdLabelFunction = vec(nbVertices);
  for (int partId = 0; partId < elementsBD().size(); partId++){ /* loop on parts */
    for (vector<int> simplex : elementsBD()[partId]){ /* loop on boundary elements */
      for (int i = 0; i < simplex.size(); i++){
        vecSetInsert(bdLabelFunction, simplex[i], bdLabels()[partId]);
      }
    }
    VecAssemblyBegin(bdLabelFunction);
    VecAssemblyEnd(bdLabelFunction);
  }
  return bdLabelFunction;
}
