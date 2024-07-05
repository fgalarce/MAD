/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2021,
    
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

#ifndef MAD_MEASURES
#define MAD_MEASURES

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
#include <innerProduct.hpp>
//#include <boundaryConditions.hpp>
#include <boundaries.hpp>
#include <linearAlgebra.hpp>

using namespace std;

class Measures {

  public:
    Measures(){}
    ~Measures(){}

    void initialize(Parameters parameters, Geometry & geometry, LinearAlgebra & la, IO & io);
    void finalize();

    inline const vector<Vec> & rieszRepresenters() const {
      return m_rieszRepresenters;}

    inline const Vec & rieszRepresenters(int i) const {
      return m_rieszRepresenters[i];}

    inline const int & nbMeasures() const {
      return m_nbMeasures;}

    inline const int & nbVoxels() const {
      return m_nbVoxels;}

    inline const double & std_dev() const {
      return m_std_dev;}

    vector<vector<double>> coordinates() {
      return m_coordinates;}

    inline const Vec & measures() const {
      return m_l;}

    inline const Vec & measuresOnMesh() const {
      return m_measuresOnMesh;}

    inline const vector<double> & rrNorms() const {
      return m_rrNorms;}

    void addNoise(double (*noise) (double,Parameters), int variableLabel = 0);
    void computeSynthetic(double time);
    void computeSynthetic(int iteration);
    void computeNoiseIntensity();
    void computeMeasuresMax();

    void computeRieszRepresentersH1();
    void computeRieszRepresentersL2();

    vector<vector<vector<double>>> matrixImage();

//    double noise_intensity;

  private:

    Parameters par;

    /* The object is designed to be initialized with a Linear Algebra */
    void initialize(Geometry geometry, InnerProduct & ip);

    /* number of sample volumes */
    int m_nbMeasures;
    int m_nbDofs;

    int m_nbVertices;

    /* frame of reference for the adta*/
    vector<vector<double>> m_imageFOR; 
    vector<vector<vector<int>>> elements;
    vector<vector<double>> m_coordinates;
    double * m_CFI; /* one-dimensional color flow image */
    vector<vector<double>> m_VFI; /* two-dimensional vector flow image */
    Vec m_l; /* CFI or VFI in petsc format. VFI is ordered as follows: v1x v2x ... vmx v1y .... vm */
    int m_nbVoxels;
    int m_nbComponentsPerVoxel;
    vector<double> m_measureMaxValue;
  
    PetscInt * indexesCFI;
    PetscInt * indexesVFI;

    vector<int> id_nz_rr;
    vector<Vec> m_rieszRepresenters;
    vector<Vec> m_rieszRepresentersL2;
    vector<double> m_rrNorms;
    vector<double> m_measures;
    Vec m_measuresOnMesh;
    Geometry geo;

    IO io_measures; // deprecated
    IO io;

    Vec u_bar_measured;

    /* Hilbert */
    InnerProduct m_ip;
    LinearAlgebra m_la;

    PetscErrorCode code;
    int m_world_rank;
    double m_std_dev;

};  

#endif
