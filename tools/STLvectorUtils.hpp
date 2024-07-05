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

#ifndef STL_VEC_UTILS 
#define STL_VEC_UTILS 

#include<math.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<iomanip>
#include<algorithm>
#include<string>
#include<typeinfo>
#include<assert.h>
#include<tools.hpp>

using namespace std;

typedef vector<double> VEC;
typedef vector<vector<double>> MAT;

MAT stl_mat(int m, int n);
int count_zeros(vector<double> u);

vector<vector<double>> identity(int dimension);

double average(vector<double> const& v);

/* return projection of u into v */
VEC proj(VEC u, VEC v, string innerProd="l2");

/* return dot product of u and v */
double dot(const vector<double> & u, const vector<double> & v);
double dot(const vector<vector<double>> & u, const vector<vector<double>> & v);

/* norm of u given an inner product */
double norm(const vector<double> & u, string innerProd="l2");
void vecNormalize(VEC u);

/* cross product between u and v */
vector<double> cross(const vector<double> & u, const vector<double> & v);

/* return largest element in u */
template <typename T>
T max(vector<T> u);

/* return sum of the elements in u */
double sum(VEC u);

/* return a vector VEC such that V[i] = abs(u[i]) */
VEC absSTL(VEC u);

double determinant(const vector<vector<double>> & m);
double determinant22(const vector<vector<double>> & m);
double determinant33(const vector<vector<double>> & m);
double determinant44(const vector<vector<double>> & m);

double trace(const vector<vector<double>> & m);

/* inverse matrix */
vector<vector<double>> invert(const vector<vector<double>> & m);

/* matrix transpose */ 
vector<vector<double>> transpose(const vector<vector<double>> & m);
vector<vector<int>> transpose(const vector<vector<int>> & m);
vector<vector<double>> transpose(const vector<double> & m);

/* Functions to import/export data */
vector<vector<double>> importdata(string fileName, int skiprows = 0);
vector<double> importdata1D(string fileName, string option = "talk");
//vector<double> importdata1d(string fileName);
//void  importdata1d(string fileName, vector<int> vecToImport);

void exportData(string filePath, vector<vector<double>> dataToExport, string delimiter=" ", int precision=6);
void exportData(string filePath, vector<vector<int>> dataToExport, string delimiter=" ");
void exportData(string filePath, vector<double> dataToExport, string delimiter = "\n", int precision=6);
void exportData(string filePath, vector<int> dataToExport);
void exportData(vector<double> dataToExport, string filePath);

/* 
    Get a specific line from an ASCII file (return error code)
*/
unsigned int getLine(string filePath, int nLine, string * theLine);

/* 
    Delete row of a double std matrix = vector<vector <double> >
*/
void deleteRow( int idx, vector<vector<double>>  vec);


/* vertcat matlab-like function for STL vectors. TODO: the same for matrices */
vector<double> vertcat( vector<double> u, vector<double> v );


/* return indices of repeated elements in u */
vector<int> findRepeatedIndices(vector<double> u);

/* erase elements of an STL vector given a vector of increasing indices */
vector<double> eraseElements(vector<double> u, vector<int> indices);

/* max and min of a STL vector in Matlab style */
double max(vector<double> u);
double min(vector<double> u);

/* find function for stl vectors. Coding:       1: =     2: <       3: >      */
vector<int> findSomething(vector<double> u, double whatToLookFor, int comparison_operator = 1);
int findSomething(const vector<int> u, int whatToLookFor);
template<class InputIterator, class T>
  InputIterator find (InputIterator first, InputIterator last, const T& val);

/* get sub-vector of a stl vector given a vector of incremental indices*/
vector<double> getSubVector(vector<double> u, vector<int> indices);

/* get some row or column of a STL matrix */
vector<double> getRow(vector<vector<double>> u, int rowId);
vector<double> getCol(vector<vector<double>> u, int colId);

/* Returns interpolated value at x from parallel arrays ( xData, yData )
   Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
   boolean argument extrapolate determines behaviour beyond ends of array (if needed)	*/
double interpolate( vector<double> &xData, vector<double> &yData, double x, bool extrapolate = false);

/* Matlab like linear interpolation function.
	INPUT:  sampling points stl vectors x y and query points xq
	OUTPUT: interpolated values over query points
*/
vector<double> interp1(vector<double> &x, vector<double> &y, vector<double> &xq);

/* double vector to string vector */
vector<string> num2Str(vector<double> u);

/* double matrix to string matrix */
vector<vector<string>> num2Str(vector<vector<double>> m);

/* string vector to double vector */
vector<double> str2Num(vector<string> u);

/* string matrix to double matrix */
vector<vector<double>> str2Num(vector<vector<string>> m);

/* horizontal concatenation of vectors matlab-style */
vector<vector<double>> horzcat(vector<double> u, vector<double> v, vector<double> w);

/* write binary file from double vector */
void writeBinaryVec(vector<double> u, string fileName);

/* write binary file from double matrix */
void writeBinaryMat(vector<vector<double>> m, string fileName);

/* read binary file to double STL matrix */
vector<vector<double>> readBinaryMat(string fileName);

/* read binary file to double STL vector */
vector<double> readBinaryVec(string fileName);

/* overload operators */
template <typename T>
ostream& operator<<(ostream& out, vector<T> u){
  for (int i=0; i<u.size(); i++){
    out << scientific << u[i] << " ";
    if (((i+1) % 6) == 0){
      out << endl;
    }
  }
  return out;
}

ostream& operator<<(ostream& out, MAT m);
ostream& operator<<(ostream& out, vector<vector<int>> m);
VEC operator/(VEC u, double alpha);

vector<vector<double>> operator/(vector<vector<double>> m, double alpha);

/* ------ Templates for vector operations ------ */
template<typename T>
vector<T> operator+(vector<T> u, vector<T> v){
  assert(u.size() == v.size());
  vector<T> sumElemWise;
  sumElemWise.resize(u.size());
  for (int i = 0; i < u.size(); i++)
    sumElemWise[i] = u[i] + v[i];
  return sumElemWise;
}

template<typename T>
vector<T> operator-(vector<T> u, vector<T> v){
  assert(u.size() == v.size());
  vector<T> sumElemWise;
  sumElemWise.resize(u.size());
  for (int i = 0; i < u.size(); i++)
    sumElemWise[i] = u[i] - v[i];
  return sumElemWise;
}

template<typename T>
vector<T> operator*(T alpha, vector<T> u){
  for (int i=0; i < u.size(); i++){
    u[i] = u[i]*alpha;
  }
  return u;
}

template<typename T>
vector<T> operator*(vector<T> u, T alpha){
  return alpha*u;
}

/* ------ Templates for matrix operations ------ */

template <typename T>
vector<vector<T>> operator+(vector<vector<T>> m, T alpha){
  vector<vector<T>> m_plus_alpha;
  m_plus_alpha.resize(m.size());
  for (int i = 0; i < m.size(); i++){
    m_plus_alpha[i].resize(m[i].size());
    fill(m_plus_alpha[i].begin(), m_plus_alpha[i].end(), 0.0);
    for (int j = 0; j < m[i].size(); j++){
      m_plus_alpha[i][j] = m[i][j] + alpha;
    }
  }
  return m_plus_alpha;
}

template <typename T>
vector<vector<T>> operator+(vector<vector<T>> m1, vector<vector<T>> m2){
  vector<vector<T>> m1_plus_m2;
  m1_plus_m2.resize(m1.size());

  assert(m1.size() == m2.size());
  assert(m1[0].size() == m2[0].size());

  for (int i = 0; i < m1.size(); i++){
    m1_plus_m2[i].resize(m1[i].size());

    for (int j = 0; j < m1[i].size(); j++){
      m1_plus_m2[i][j] = m1[i][j] + m2[i][j];
    }
  }
  return m1_plus_m2;
}

template <typename T>
vector<vector<T>> operator-(vector<vector<T>> m, T alpha){
  vector<vector<T>> m_minus_alpha;
  m_minus_alpha.resize(m.size());
  for (int i = 0; i < m.size(); i++){
    m_minus_alpha[i].resize(m[i].size());
    fill(m_minus_alpha[i].begin(), m_minus_alpha[i].end(), 0.0);
    for (int j = 0; j < m[i].size(); j++){
      m_minus_alpha[i][j] = m[i][j] - alpha;
    }
  }
  return m_minus_alpha;
}

template<typename T>
vector<vector<T>> operator*(vector<vector<T>> m, T alpha){
  vector<vector<T>> m_times_alpha;
  m_times_alpha.resize(m.size());
  for (int i = 0; i < m.size(); i++){
    m_times_alpha[i].resize(m[i].size());
    fill(m_times_alpha[i].begin(), m_times_alpha[i].end(), 0.0);
    for (int j = 0; j < m[i].size(); j++){
      m_times_alpha[i][j] = m[i][j]*alpha;
    }
  }
  return m_times_alpha;
}

template<typename T>
vector<vector<T>> operator*(T alpha, vector<vector<T>> m){
  return m*alpha; /* forward due to commutative operator */
}

template<typename T>
vector<vector<T>> operator*(vector<vector<T>> m1, vector<vector<T>> m2){

  int nbRows = m1.size();
  int nbCols = m2[0].size();

  vector<vector<T>> m1_times_m2(nbRows);

  for (int i = 0; i < nbRows; i++){
    m1_times_m2[i].resize(nbCols);
    for (int j = 0; j < nbCols; j++){
      m1_times_m2[i][j] = dot(getRow(m1, i), getCol(m2, j));
    }
  }

  return m1_times_m2;
}

template<typename T>
vector<T> operator*(vector<vector<T>> m, vector<T> u){
  assert(m[0].size() == u.size());
  int nbRows = m.size();
  vector<T> m_times_u(nbRows);
  for (int i = 0; i < nbRows; i++){
    m_times_u[i] = dot(m[i], u);
  }
  return m_times_u;
}

//template<typename T>
//vector<T> operator*(vector<T> u, vector<vector<T>> m){
//  return m*u;
//}

inline vector<int> range(int n){
  vector<int> np_range(n);
  for (int i = 0; i < n; i++){
    np_range[i] = i;
  }
  return np_range;
}

inline vector<int> rango(int n1, int n2){
  vector<int> index(n2 - n1 + 1);
  for (int i = 0; i < n2 - n1 + 1; i++){
    index[i] = n1 + i;
  }
  return index;
}

inline vector<double> rango(double x_start, double x_end, double step){
  vector<double> np_array((x_end - x_start + step)/step);
  for (int i = 0; i < np_array.size(); i++){
    np_array[i] = x_start + i*step; 
  }
  return np_array;
}

inline vector<int> range(int n1, int n2){
  vector<int> index(n2 - n1);
  for (int i = 0; i < n2 - n1; i++){
    index[i] = n1 + i;
  }
  return index;
}

inline vector<double> range(double x_start, double x_end, double step){
  vector<double> np_array((x_end - x_start)/step);
  for (int i = 0; i < np_array.size(); i++){
    np_array[i] = x_start + i*step; 
  }
  return np_array;
}

void erase(vector<double> & u, vector<int> indices);

void erase(vector<vector<double>> & u, vector<int> indices);

vector<double> arcsin(vector<double> u);
vector<double> arccos(vector<double> u);

#endif
