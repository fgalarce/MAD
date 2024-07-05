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

#ifndef ULTRA_4D_FLOW_MODEL
#define ULTRA_4D_FLOW_MODEL

#include <iostream>
#include <math.h>
#include <vector>
#include <slepc.h>
#include <STLvectorUtils.hpp>
#include <petscWrapper.hpp>
#include <tools.hpp>
#include <parameters.hpp>
#include <geometry.hpp>
#include <linearAlgebra.hpp>

using namespace std;

class Model {

  public:
    Model(){}
    ~Model(){}

    void initialize(Parameters parameters, size_t nbVertices);
    void initialize(Parameters parameters, Geometry & geometry, LinearAlgebra & la, IO & inputOutput);
    void finalize();

    const vector<Vec> & basis(const double time); /* Constant reference returned from a constant function. A constant function is a function of a class with no right to make change in the class objects */

    const Vec average();

    inline const bool & updateBasis() const {
      return m_updateBasis; }
    
    inline const int & indexCurrentBasis() const {
      return m_indexCurrentBasis;}

    inline const int & coord_HR() const {
      return m_coord_i;}

    inline const vector<int> & nbModesMapping() const {
      return m_nbModesMapping;}

    inline const vector<double> & delta() const {
      return m_delta[m_indexCurrentBasis];}

  private:

    LinearAlgebra m_la;
    IO io;

    void initialize_PBDW_pice_wise();
    void initialize_PBDW_linear();

    vector<vector<Vec>> m_basis;

    Parameters par;

    int  m_nbVertices;

    int m_coord_i; /* HR coordinate in data-set */

    int m_world_rank;

    vector<Vec> m_uBar; /* average of snapshots at each time window */

    int m_indexCurrentBasis = 0;  /* index of time-window of current reconstruction */
    int m_indexPreviousBasis = -1;
    bool m_updateBasis = true;

    Geometry geo;

    vector<int> m_nbModesMapping;
    vector<PetscInt> m_indexes;
    vector<PetscInt> m_indexes_vel;
    vector<PetscInt> m_indexes_press;
    vector<PetscInt> m_indexes_press_mix;

    /* constrained least-squares to deal with noise */
    vector<vector<double>> m_delta;

};  

#endif
