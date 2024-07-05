/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2024,
    
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

#ifndef MAD_GEOMETRY
#define MAD_GEOMETRY

#include <iostream>
#include <math.h>
#include <iterator>
#include <vector>
#include <slepc.h>
#include <STLvectorUtils.hpp>
#include <petscWrapper.hpp>
#include <tools.hpp>
#include <parameters.hpp>
#include <string>
#include <io.hpp>
#include <masterElement.hpp>

using namespace std;

class Geometry {

  public:
    Geometry(){}
    ~Geometry(){}
    void initialize(Parameters parameters, const IO & io);
    void initialize(const IO & io);
    void initialize(const IO & io, vector<double> T); /* To be used wisely :)*/

    int nearest_neighbourd(vector<double> point); 

    /* check if point is inside a box */
    bool isInside(const vector<double> & point, const vector<double> & voxelCenter, 
                  const vector<double> & beamDir, const vector<double> & transDir, const vector<double> & zDir, 
                  double voxelSizeX, double voxelSizeY, double voxelSizeZ);
    bool isInside(const vector<double> & point, const vector<double> & voxelCenter, 
                  const vector<double> & beamDir, const vector<double> & transDir, const vector<double> & zDir);
    bool isInside(const vector<double> & point, const vector<double> & voxelCenter, const vector<vector<double>> & voxelBasis,
                  double voxelSizeX, double voxelSizeY, double voxelSizeZ); 
    bool isInside(const vector<double> & point, const vector<double> & voxelCenter, const vector<vector<double>> & voxelBasis, 
                  const vector<double> & voxelSize); 

    inline const vector<vector<double>> & coordinates() const {
      return m_coordinates;}

    inline const vector<int> & coordinatesLabels() const {
      return m_coordinatesLabels;}

    inline const vector<double> & coordinates(int i) const {
      return m_coordinates[i];}

    inline const vector<int> & edgeLabels() const {
      return m_edgeLabels;}

    inline const int & edgeLabels(int i) const {
      return m_edgeLabels[i];}

    inline const vector<int> & triangleLabels() const {
      return m_triangleLabels;}

    inline const int & triangleLabels(int i) const {
      return m_triangleLabels[i];}

    inline const int & dimension() const {
      return m_dimension;}

    inline const vector<int> & tetrahedronLabels() const {
      return m_tetrahedronLabels;}

    inline const vector<vector<vector<int>>> & tetrahedron() const{  
      return m_tetrahedron;}

    inline const vector<vector<vector<int>>> & triangles() const {  
      return m_triangles;}

    inline const vector<vector<vector<int>>> & edges() const { 
      return m_edges;}

    inline const vector<vector<int>> & boundaryNodes() const { 
      return m_boundaryNodes;}

    inline const vector<int> & boundaryNodes(int label) const { 
      if (m_dimension == 2) {
        return m_boundaryNodes[findSomething(m_edgeLabels, label)];
      } else {
        return m_boundaryNodes[findSomething(m_triangleLabels, label)];
      }
    }

    inline const vector<int> & bdLabels() const { 
      if (m_dimension == 2) {
        return m_edgeLabels;
      } else {
        return m_triangleLabels;
      }
    }

    inline const vector<vector<int>> & tetrahedron(int label) const{ 
      return m_tetrahedron[findSomething(m_tetrahedronLabels, label)];}

    inline const vector<vector<int>> & quad(int label) const{ 
      return m_quad[findSomething(m_quadLabels, label)];}

    inline const vector<vector<int>> & triangles(int label) const {  
      return m_triangles[findSomething(m_triangleLabels, label)];}

    inline const vector<vector<int>> & edges(int label) const { 
      return m_edges[findSomething(m_edgeLabels, label)];}

    inline const double & meshVolume() {
      return m_meshVolume;
    }

    inline const int & nbTetrahedra() const {
      return m_nbTetrahedra;}

    inline const int & nbTriangles() const {
      return m_nbTria;}

    inline const int & nbNodes() const {
      return m_nbNodes;}

    inline const int & nbElementsBoundary() const {
      return m_nbElementsBoundary;}

    inline const int & nbBoundaryNodes() const {
      return m_nbBoundaryNodes;}

    int nbVertices;

    vector<vector<int>> getElementNeighbourhood(int pointId);
    vector<vector<int>> getTetraNeighbourhood(int pointId);
    vector<int> star(int pointId);
    vector<vector<int>> getTriaNeighbourds(int pointId);
    vector<vector<int>> getTriaNeighbourdsDEPRECATED(int pointId, int i);
    vector<int> getPointStarBD(int pointId, int i);
    vector<vector<int>> getPointStarBD(int pointId);
    vector<int> getTetraFromTriangle(vector<int> triangle); /* to be deprecated, as every dimension dependent function */
    vector<int> getElementFromElementBD(vector<int> simplexBD);
    vector<double> getNodalValues(Vec u0_sequential, vector<int> simplex, int nbDofsPerNode);

    void computeMeshVolume();
    vector<double> getCenter(int labelBD);
    double getRadius(int labelBD, vector<double> center);
    double volume();
    double computeBoundaryArea();
    double computeBoundaryArea(int labelBD);
    void print_geo();
    Vec boundaryLabelFunction();

    inline const vector<vector<int>> & elements(int label) const{ 
      return m_elements[findSomething(m_elementLabels, label)];}

    inline const vector<vector<int>> & elementsBD(int label) const{ 
      return m_elementsBD[findSomething(m_elementLabelsBD, label)];}

    inline const int & elementLabels(int i) const {
      return m_elementLabels[i];}

    inline const int & elementLabelsBD(int i) const {
      return m_elementLabelsBD[i];}

    inline const vector<int> & elementLabelsBD() const {
      return m_elementLabelsBD;}

    inline const vector<vector<vector<int>>> & elements() const{  
      return m_elements;}

    inline const vector<vector<vector<int>>> & elementsBD() const{  
      return m_elementsBD;}

    inline const Mat & massMatrix() const{
      return m_M;}

  private:
    Parameters par;

    vector<vector<vector<int>>> m_elements;
    vector<vector<vector<int>>> m_elementsBD;
    vector<int> m_elementLabels;
    vector<int> m_elementLabelsBD;

    void assembleMass();
    void deployGeo(const IO & io);
    void deployGeo1D(const IO & io);
    void deployGeo2D(const IO & io);
    void deployGeo3D(const IO & io);
    vector<vector<double>> m_coordinates;
    vector<vector<vector<int>>> m_tetrahedron; /* part (subdomain label) -> tetra -> points ids */
    vector<vector<vector<int>>> m_triangles; /* part (subdomain label) -> tria -> points ids */
    vector<vector<vector<int>>> m_edges; /* part (subdomain label) -> edge -> points ids */
    vector<vector<vector<int>>> m_quad; /* part (subdomain label) -> tria -> points ids */
    vector<vector<int>> m_boundaryNodes; /* part (boundary label) -> point ids */

    vector<int> m_coordinatesLabels;
    vector<int> m_edgeLabels;
    vector<int> m_quadLabels;
    vector<int> m_triangleLabels;
    vector<int> m_tetrahedronLabels;

    double m_meshVolume;

    bool m_mesh_volume_is_calculated = false;

    int m_dimension;
    int m_nbTria;
    int m_nbTetrahedra;
    int m_nbNodes;
    int m_nbElementsBoundary;
    int m_nbBoundaryNodes;

    int m_world_rank;
    MasterElement fe;
    Mat m_M;

    int m_verbose = 4;
    PetscErrorCode code;
};  

#endif
