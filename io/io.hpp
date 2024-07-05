/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2020,
    
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

#ifndef MAD_IO
#define MAD_IO

/* 
 * Handle the read and the writing of meshes, solutions and time stepping.
 * Base format is that of Ensight. There is support for INRIA .mesh input
 *
 * Felipe Galarce 2020
 */ 

#include <iostream>
#include <fstream>
#include <parameters.hpp>
#include <vector>
#include <STLvectorUtils.hpp>
#include <tools.hpp>
#include <iomanip> /* setprecision(int) */
#include <slepc.h>
#include <typeinfo>
#include <petscWrapper.hpp>

using namespace std;

class IO {

  public:

    IO(){}

    void initialize(string geometryData, string dirResults);
    void initialize(Parameters parameters, string dirSurname = "", string jobSurname = "");
    void initialize(string geometryData, string dirResults, Parameters parameters);

    void loadVector(Vec u, string filename);
    void loadVector(Vec u, vector<string> filename);
    void loadVector(Vec u, string filename0, string filename1);

    /* Write state per node */
    void writeState(double * u, string variableName, int vecSize, double time = -1.0);
    void writeState(vector<double> u, string variableName, double time = -1.0); /* wrapper */
    void writeState(Vec u, double time = -1.0, const vector<string>& varName= vector<string>()); /* wrapper */
    void writeState(Vec u, string variableName, double time = -1.0); /* wrapper */

    /* Write state per element */
    void writeStateP0(Vec u, string variableName, double time = -1.0); /* wrapper */
    void writeStateP0(vector<double> u, string variableName, double time = -1.0); /* wrapper */
    void writeStateP0(double * u, string variableName, int vecSize, double time = -1.0);

    /* VTK */
    void writeVTKgrid(vector<double> data, int nx, int ny, int nz, string name, 
                      double aspect_ratio_x = 1.0, double aspect_ratio_y = 1.0, double aspect_ratio_z = 1.0);
    /* Acces functions */
    inline const vector<vector<double>> & coordinates() const {
      return m_coordinates;}

    inline const vector<int> & coordinatesLabels() const {
      return m_coordinatesLabels;}

    inline const vector<int> & edgeLabels() const {
      return m_edgeLabels;}

    inline const vector<int> & triangleLabels() const {
      return m_triangleLabels;}

    inline const vector<int> & quadLabels() const {
      return m_quadLabels;}

    inline const vector<int> & tetrahedronLabels() const {
      return m_tetrahedronLabels;}

    inline const vector<vector<vector<int>>> & tetrahedron() const {
      return m_tetrahedron;}

    inline const vector<vector<vector<int>>> & triangles() const {
      return m_triangles;}

    inline const vector<vector<vector<int>>> & edges() const {
      return m_edges;}

    inline const int & nbVertices() const{
      return m_nbVertices;}

    inline const int & nbTriangles() const{
      return m_nbTriangles;}

    inline const int & dimension() const{
      return m_dimension;}

    inline const int & nbElementsBoundary() const{
      return m_nbElementsBoundary;}

    Vec P1bP1_to_P1P1(Vec u);
  private:

    Parameters par;

    void readGeo();
    void readGeoEnsight();
    void readGeoEnsightOld();
    void readGeoInria();
    void readGeoInria1D();
    void readGeoInria2D();
    void readGeoInria3D();
    void writeGeo();

    vector<double> m_time;

    string m_jobName;
    string m_outputDir;
    string m_geometryName;

    vector<string> m_variables;
    vector<string> m_variablesType;
    vector<string> m_variablesTypeNodeElement;

    /* geometry */
    string m_geometryData;
    string m_nameGeometry;
    int m_nbVertices;
    int m_nbEdges;
    int m_nbTriangles;
    int m_nbElementsBoundary;
    int m_nbQuad;
    int m_nbTetrahedron;
    vector<vector<double>> m_coordinates;
    vector<vector<vector<int>>> m_edges; /* part (subdomain label) -> triangle -> points ids */
    vector<vector<vector<int>>> m_triangles; /* part (subdomain label) -> triangle -> points ids */
    vector<vector<vector<int>>> m_quads; /* part (subdomain label) -> tetra -> points ids */
    vector<vector<vector<int>>> m_tetrahedron; /* part (subdomain label) -> tetra -> points ids */

    vector<int> m_coordinatesLabels;
    vector<int> m_edgeLabels;
    vector<int> m_triangleLabels;
    vector<int> m_quadLabels;
    vector<int> m_tetrahedronLabels;
  
    void writeCase();

    void updateCase(string variable, string typeVariable, string perNodeOrPerElement, double time);

    string getElementType(int partId);

    int m_world_rank;

    int m_dimension;

    PetscErrorCode code;
};  

#endif
