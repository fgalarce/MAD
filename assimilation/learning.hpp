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

#ifndef ULTRA_4D_FLOW_LEARNING
#define ULTRA_4D_FLOW_LEARNING

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <random>
#include <chrono>
#include <slepc.h>
#include <STLvectorUtils.hpp>

#include <petscWrapper.hpp>
#include <tools.hpp>
#include <parameters.hpp>
#include <geometry.hpp>
#include <io.hpp>
//#include <masterTetrahedron.hpp>

using namespace std;

class Learning {

  public:
    Learning(){}
    ~Learning(){}

    void initialize(Parameters parameters, Geometry geometry);
    void initialize(Parameters parameters);

    inline const vector<double> & voxelized_geometry() const {
      return m_voxelized_geometry;}

    void voxelize();
    void write_geo_unstructured_grid();
    Mat & MDS(Mat D);

    /* Least squares */
    void assemble_normal_equations();
    vector<double> compute_weights(int idRow);
    inline const double & residual() const { 
      return m_residual;}
    void assembleG();

  private:

    Parameters par;
    PetscErrorCode code;

    int m_nbDofs, m_nbVertices;

    void get_clouds_per_voxel();
    void voxelize_points();
    void voxelize_volume();

    int m_nbVoxels;
    vector<vector<double>> m_voxels;
    vector<vector<int>> m_clouds;
    Geometry m_geo;

    int m_world_rank;

    vector<double> m_voxelized_geometry;
    vector<vector<double>> dissimilarityMatrix(vector<vector<double>> Y);

//    int m_nbSamples;

    /* Least squares */
    Mat m_G, m_pinvG;
    double m_residual;
    Vec m_w;

};  

#endif
