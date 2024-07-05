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

#ifndef MAD_BOUNDARIES
#define MAD_BOUNDARIES

#include <iostream>
#include <math.h>
#include <vector>
#include <slepc.h>
#include <STLvectorUtils.hpp>
#include <petscWrapper.hpp>
#include <tools.hpp>
#include <parameters.hpp>
#include <masterElementBD.hpp>
#include <io.hpp>
#include <geometry.hpp>
#include <interpolate.hpp>

using namespace std;

class Boundary{

  public:
    Boundary(){}
    ~Boundary(){}

    void initialize(Parameters parameters, const Geometry & geometry);

    inline const Vec normals() const { /* nbDofs size */
      return m_normals;}

    inline const Vec normals_only() const { /* m_dimension * geo.nbVertices size*/
      return m_normals_only;}

    inline const vector<vector<double>> normalsP0comp() const {
      return m_normals_P0comp;}

    inline const vector<double> normalsP0() const {
      return m_normals_P0;}

    inline const vector<double> paraboloid_center(int bdLabel){
      return m_paraboloid_center[findSomething(geo.elementLabelsBD(), bdLabel)];}

    inline const double paraboloid_radius(int bdLabel){
      return m_paraboloid_radius[findSomething(geo.elementLabelsBD(), bdLabel)];}

    inline const Mat NeumannMass(int bdLabel){
      return m_NeumannMass[findSomething(geo.elementLabelsBD(), bdLabel)];}

    inline const Mat scalarMass(int bdLabel){
      return m_scalarMass[findSomething(geo.elementLabelsBD(), bdLabel)];}

    inline const Mat mass(){
      return m_M;
    }

    /* Dirichlet bc */
    void Dirichlet                (int bdLabel, Vec u, int variableLabel = 0); 
    void Dirichlet                (int bdLabel, vector<double> (*bc_func) (vector<double>,double,Parameters), int variableLabel = 0);
    void Dirichlet_2D_or_3D       (int bdLabel, vector<double> (*bc_func) (vector<double>,double,Parameters), int variableLabel = 0);
    void DirichletComp(int compBC, int bdLabel, vector<double> (*bc_func) (vector<double>,double,Parameters), int variableLabel = 0);
    void DirichletParaboloid      (int bdLabel, vector<double> (*bc_func) (vector<double>,double,Parameters), int variableLabel = 0);
    void DirichletNormalParaboloid(int bdLabel,        double  (*bc_func) (vector<double>,double,Parameters), int variableLabel = 0);
    void DirichletNormalParaboloid(int bdLabel,        double  (*bc_func) (vector<double>,double,Parameters), vector<int> variableLabel);
    void reflection(Vec u, int bdLabel);

    /* Neumann bc */
    void Neumann                (int bdLabel, vector<double> (*bc_func) (vector<double>,double,Parameters), int variableLabel = 0);
    void Neumann_2D_or_3D       (int bdLabel, vector<double> (*bc_func) (vector<double>,double,Parameters), int variableLabel = 0);
    void NeumannNormal          (int bdLabel,        double  (*bc_func) (vector<double>,double,Parameters), int variableLabel = 0);
    void NeumannNormalParaboloid(int bdLabel,        double  (*bc_func) (vector<double>,double,Parameters), int variableLabel = 0);
    void NeumannNormalConstant  (int bdLabel,        double  bc_value);
    void NeumannNormalBackflow  (int bdLabel, Vec stab_force); // HAML added

    /* methods to impose bc */
    void block(Mat & A, string method="brut");
    void block(Vec & rhs, string method="brut");
    void block(Mat & A, Vec & rhs, string method="brut");

    /* Calculator */
    double computeSize(int bdLabel);

    /* Stabilization */
    void backflowNeumann(Vec u, int bdLabel); // HAML added
    void backflow(Vec u);
    void backflowStrong(Vec u);
    void backflowStrong(Vec u, int idBD);

    /* Mass */
    MasterElementBD feBD;
    void assembleBoundaryMass(string mode = "contiguous", int nbVar = -1);
    void assembleBoundaryNormals(string mode = "contiguous");
    void assembleBoundaryMassScalar();

    double flow(Vec u, int labelBD);
    double integral(Vec p, int labelBD);

    Parameters par;

    PetscErrorCode code;
    int m_verbose = 4;

    int m_world_rank;

    Mat m_M;

    vector<Mat> m_scalarMass;
    vector<Mat> m_NeumannMass;

    void time(double t);

  private:
    double m_time = 0.0;
    int m_iteration = 0;
    vector<vector<double>> m_inletU0;
    Geometry geo;

    int nbDofs;
    int m_dimension;
    vector<vector<int>> m_idDofsBC; /* part --> id dofs*/
    vector<vector<double>> m_bc;    /* part --> value vc at id m_idDofsBC*/
    Vec m_NeumannBC_weak;
    void reset();

    /* Geometrical information */
    Vec m_normals; /* with nbDofs size */
    Vec m_normals_only; /* with m_dimension*geo.nbVertices size */
    Vec m_normals_sequential;
    vector<double> m_stl_normals_sequential;
    vector<vector<double>> m_paraboloid_center;
    vector<double> m_paraboloid_radius;
    vector<vector<double>> m_normals_P0comp;
    vector<double> m_normals_P0;
    vector<int> m_corners;

    void assembleNormals2D();
    void computeCorners();
      
    bool shut_up = true;
};  

#endif
