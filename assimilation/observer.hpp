/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
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

#ifndef MAD_OBSERVER
#define MAD_OBSERVER

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
#include <innerProduct.hpp>
#include <boundaries.hpp>
#include <linearAlgebra.hpp>
#include <calculus.hpp>

using namespace std;

class Observer {

  public:
    Observer(){}
    ~Observer(){}

    void initialize(Parameters parameters, Geometry & geometry, InnerProduct & innerProduct, IO & inputOutput, Calculus & calc, LinearAlgebra & linearAlgebra, Boundary & boundary);
    void finalize();

    Vec synthetic(int iteration);
  
    void buildObserver();
    void buildMeasures(int iteration);
    void exportObserver();
    void loadData();
    vector<Vec> observations(int iteration);
    vector<Vec> observations(int iteration, string dirField, string sufix);
    Vec loadExpected(string expec_folder);
    void printRieszRepresenters();

    inline const Vec basis(int i) const {
      return m_rieszRepresenters[i];}

    inline const Vec basis_l2(int i) const {
//      Vec basis_l2 = vec(m_rieszRepresenters_l2[i]);
//      VecScale(basis_l2, m_rrNorms_l2[i]); 
      return m_rieszRepresenters_l2[i];}
    
    inline const double l2_rrNorm(int i){
      return sqrt(ip(m_rieszRepresenters_l2[i], m_rieszRepresenters_l2[i])); 
    }

    inline const int nbMeasures() const {
      return m_nbMeasures;}

    inline const Vec measures() const {
      return m_theta;}

  private:

    Parameters par;
    PetscErrorCode code;

    int m_nbDofs;

    vector<Vec> m_rieszRepresenters;
    vector<Vec> m_rieszRepresenters_l2;
    vector<double> m_rrNorms;
    vector<double> m_rrNorms_l2;

    vector<vector<double>> m_voxelCenters;
    vector<vector<double>> m_voxelBasis;
    vector<vector<double>> m_measures;
    int m_nbVoxels;

    int m_world_rank;

    Geometry geo;
    InnerProduct ip;
    LinearAlgebra la;
    IO io;
    Calculus calculus;
    MasterElement fe; 
    Boundary bd;

    int m_nbMeasures;
    Vec m_omega; /* measures coordinates in V */
    Vec m_theta; /* measures coordinates in W_m */
  
    vector<int> zeroVoxels;
    vector<int> nzVoxels;
    vector<vector<double>> data;
    vector<vector<int>> m_pointsInside;

    void build_l2_Observer(int idVar, int offset);
    void buildObserverFromIP(int idVar, int offset);
    KSP ksp;
};  

vector<double> noslip_bc(vector<double> x, double t, Parameters par){
  vector<double> bc(1, 0.0);
  return bc;
}

#endif
