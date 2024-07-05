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

#ifndef MAD_petscwrapper
#define MAD_petscwrapper

#include <iostream>
#include <math.h>
#include <vector>
#include <slepc.h>
#include <STLvectorUtils.hpp>
#include <tools.hpp>
#include <fstream>
#include <parameters.hpp>

void vec(Vec & u, int sizeVec);
Vec vec(int sizeVec);
Vec vec(int sizeVec, string path_vec);
vector<Vec> vec(vector<Vec> u);
Vec vec(Vec u);
Vec zeros(int nbZeros);
void zeros(Vec & u, int nbZeros);

void matSet(Mat A, int i, int j, double Aij);
void matSetInsert(Mat A, int i, int j, double Aij);
void vecSet(Vec u, int i, double ui);
void vecSetInsert(Vec u, int i, double ui);

Mat eye(int n, string type = "sparse");
Vec ones(int n);
void mat(Mat & A, int nbRows, int nbCols, string type = "sparse");
Mat mat(int nbRows, int nbCols, string type = "sparse");
Mat mat(Mat A);
Mat inverse(Mat A);
Mat diag(vector<double> d);
Mat transpose(Mat A);

void saveVec(Vec u, string filename);
void saveMat(Mat A, string filename);
void loadMat(Mat M, string filename);
void loadVec(Vec u, string filename);
void loadVec(vector<double> & u, string filename);

void configureKSP(KSP & ksp, Parameters par); 
void configureKSP(KSP & ksp, Mat & A, 
              string solver = "gmres", 
              string preconditioner = "asm", 
              bool monitor = false, 
              bool use_solution_as_guess_KSP = false,
              bool reuse_preconditioner = false,
              double ksp_tolerance = 1e-5 );
void configureKSP(KSP & ksp, 
              string solver = "gmres", 
              string preconditioner = "asm", 
              bool monitor = false, 
              bool use_solution_as_guess_KSP = false, 
              bool reuse_preconditioner = false,
              double ksp_tolerance = 1e-5);
KSP configureKSP(Mat & A, string solver = "gmres", string preconditioner = "asm", bool monitor = false);

void configureEPS(Parameters par, EPS eps, Mat A);

vector<double> getSingularValues(Mat A, int nbModes);

void CHKERR(PetscErrorCode code);
PetscErrorCode MyKSPMonitor(KSP ksp,int n, double rnorm, void *dummy);
PetscErrorCode krylovMonitor(KSP ksp,int n, double rnorm, void *dummy);

double dot(Vec u, Vec v);
double norm(Vec u);
double norm(Mat m);

ostream& operator<<(ostream& out, Vec u);
ostream& operator<<(ostream& out, Mat m);
//Mat operator*(Mat A, Mat B);

double ip(Vec u, Mat M, Vec v);

void blockMatrix(Mat & A, vector<int> index);

Mat buildProjector(const vector<Vec> & basis);
vector<vector<double>> eigenValueProblem(Mat & A, Mat & B, string mode = "GNHEP");
double quadraticForm(Vec a, Mat M, Vec b);

Mat outer(Vec u);

Vec petsc(vector<double> u);
Mat petsc(vector<vector<double>> m);
vector<double> stl(Vec u);
vector<vector<double>> stl(Mat A);
Mat bestATAinverse(Mat A);
Mat bestAATinverse(Mat A);

Vec getSequential(Vec u);
//vector<Vec> split(Vec u, vector<int> nbDofsPerNode, int nbVertices);
//void split(vector<Vec> split_vector, Vec u, vector<int> nbDofsPerNode, int nbVertices);
#endif
