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

#ifndef MAD_IP
#define MAD_IP

#include <slepc.h>
#include <petscWrapper.hpp>
#include <parameters.hpp>
#include <geometry.hpp>
#include <boundaries.hpp>
#include <masterElement.hpp>
#include <fem.hpp>

using namespace std;

class InnerProduct{

  public:
    InnerProduct(){}
    ~InnerProduct(){}

    inline const Mat ip_mat() const { 
      return m_ip;}

//    void initialize(Parameters parameters, const Geometry & geo);
//    void initialize(Parameters parameters, Geometry geo, const Mat & M);
    void initialize(Parameters parameters, const Geometry & geometry, const Boundary & boundary);
    void finalize();

    double operator ()(const Vec & u, const Vec & v, int varLabel);
    double operator ()(const Vec & u, const Mat & M, const Vec & v);
    double operator ()(const Vec & u, const Vec & v);

    double only_grad(const Vec & u, const Vec & v);

    /* Acces functions */
    Mat & mass() {
      return m_mass;}

    Mat & stiff() {
      return m_stiff;}

    Mat matrix;
    vector<Mat> subMatrix; /* if more than one variable */

  private:
    Parameters par;

    /* ---- old functions ---- */
    Mat m_mass, m_stiff;
    int m_matSize, m_nbVertices;
    vector<int> m_nbDofsPerNode;
    double loadWeightPressure();
    void assembleAmbientSpaceMatrix();
    void loadAmbientSpaceMatrix();
    Geometry m_geo;
    /* ---- old functions ---- */

    FEM fem;
    void assembleInnerProduct();

    PetscErrorCode code;
    int m_world_rank, m, n;
    int m_dimension;
    int nbDofs;

    Geometry geo;
    Boundary bd;
    vector<int> nbDofVar;
    int nbVertices;
    int nbNodesPerElement;

    MasterElement fe;    
    Mat m_ip;

    int m_verbose = 4;
};  

#endif
