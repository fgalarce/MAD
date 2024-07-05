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

#include"STLvectorUtils.hpp"

VEC absSTL(VEC u){
  for (int i=0; i<u.size(); i++){
    if (u[i] > 0)
      u[i] = u[i];
    else
      u[i] = -u[i];
  }
  return u;
}

VEC proj(VEC u, VEC v, string innerProd){
  if (innerProd == "l2")
    return ( dot(u,v)/dot(v,v) )*v; 
  else {
    errorMessage("double<double> proj(vector<double> u, vector<double> v, string innerProd)", "Inner product " + innerProd + " does not exist.\n");
    return u; // this line is here only to avoid compiler warning
  }
}

vector<vector<double>> identity(int dimension){
  vector<vector<double>> I(dimension);
  for (int i = 0; i < dimension; i++){
    I[i].resize(dimension, 0.0);
    I[i][i] = 1.0;
  }
  return I;
}

MAT stl_mat(int m, int n){
  vector<vector<double>> Matrix(m);
  for (int i = 0; i < m; i++){
    Matrix[i].resize(n);
  }
  return Matrix;
}

int count_zeros(vector<double> u){
  int nnz = 0;
  for (int i=0; i<u.size(); i++){
    if (u[i] == 0) nnz = nnz + 1;
  }
  return nnz;
}

double dot(const vector<double> & u, const vector<double> & v){
  double dotP = 0.0;
  assert(u.size() == v.size());
  for (int i = 0; i < u.size(); i++){
    dotP = dotP + u[i]*v[i];
  }
  return dotP;
}

double dot(const vector<vector<double>> & u, const vector<vector<double>> & v){
  assert(u.size() == 1 && v.size() == 1 && u[0].size() == v[0].size());
  double dotP = 0.0;
  for (int i = 0; i < u[0].size(); i++){
    dotP = dotP + u[0][i]*v[0][i];
  }
  return dotP;
}

double norm(const vector<double> & u, string innerProd){
  if (innerProd == "l2"){
    return sqrt(dot(u,u)); 
  } else { 
    errorMessage("double<double> norm(vector<double> u, string innerProd)", "Inner product " + innerProd + " not valid.\n") ;
    return 0.0;
  }
}

void vecNormalize(VEC u){
  double normVec = norm(u);
  for (unsigned int i = 0; i < u.size(); i++){
    u[i] = u[i]/normVec;
  }
}

vector<double> cross(const vector<double> & u, const vector<double> & v){
  vector<double> crossProduct = {u[1]*v[2] - u[2]*v[1], u[2]*v[0] - u[0]*v[2], u[0]*v[1] - u[1]*v[0]};
  return crossProduct;
}

template <typename T>
T max(vector<T> u){
  try {
    // this works for ubuntu
    return *(max(u.cbegin(), u.cend()));
  } catch(...) {
    // this works for MACOSX
    return *(max_element(u.cbegin(), u.cend()));
  }
}

double sum(VEC u){
  double sum = 0;
  for (auto it=u.cbegin(); it!=u.cend(); it++)
    sum=sum + *it;
  return sum;
}

ostream& operator<<(ostream& out, MAT m){
  for (unsigned int i = 0; i < m.size(); i++){
    for (auto it = m[i].cbegin(); it!= m[i].cend(); it++){
      if (*it >= 0){
        out << " " << scientific << *it << " ";
      } else {
        out << scientific << *it << " ";
      }
    }
    out << endl;
  }
  return out;
}

ostream& operator<<(ostream& out, vector<vector<int>> m){
  out << "[ ";
  for (unsigned int i = 0; i < m.size(); i++){
    for (auto it = m[i].cbegin(); it!= m[i].cend(); it++)
      out << *it << " ";
    out << " ]\n";  
  }
  return out;
}

VEC operator/(VEC u, double alpha){
  VEC result;
  result = (1/alpha)*u;
  return result;
}


VEC operator/(double alpha, VEC u){
  return (1/alpha)*u;
}


vector<vector<double>> operator/(vector<vector<double>> m, double alpha){
  vector<vector<double>> m_alpha(m.size());
  for (int i = 0; i < m.size(); i++){
    m_alpha[i].resize(m[i].size());
    for (int j = 0; j < m[i].size(); j++){
      m_alpha[i][j]=m[i][j]/alpha;
    }
  }
  return m_alpha;
}

vector<vector<double>> importdata(string fileName, int skiprows) {

//  cout << "Reading: " << fileName << endl;

  vector<vector<double>> vecToImport;
  string line;
  vector<double> row;
  ifstream theFile ((fileName).c_str());

  /* jump header */
  for(int i = 0 ; i < skiprows; i++){
    getline(theFile,line);
  }

  if (theFile.is_open()){
    while (getline(theFile,line)){ 
      while (line.length() > line.find(" ")){
        /* this if/break statements  avoid stod() exception when importing files with white spaces at the end of each row */
        if ( line.substr(line.find(" ") + 1, line.length()).length()  != 0 ){
          row.push_back(stod(line.substr(0, line.find(" "))));
          line = line.substr(line.find(" ") + 1, line.length());
        } else {
          break;
        }
      }
      row.push_back(stod(line));      
      vecToImport.push_back(row);
      row.erase(row.cbegin(), row.cend());
    }
  } else {
    errorMessage("importdata", "the file "+fileName+" cannot be opened or does not exists.\n");
  }
  return vecToImport;
}


vector<double> importdata1D(string fileName, string option){

  if (option != "shut_up"){
    cout << "Reading: " << fileName << endl;
  }

  /* One column data
    v1
    v2
    ...
    vn
  */
  vector<double> vecToImport;

  string line;

  ifstream theFile (fileName);

  if (theFile.is_open()){
    while (getline(theFile,line)){ 
      vecToImport.push_back(stod(line));
    }
  } else {
    errorMessage("importdata", " the file "+fileName+" cannot be opened or does not exists.\n");
  }
  return vecToImport;
}

//vector<double> importdata1d(string fileName){
//
//  vector<double> vecToImport;
//
//  string line;
//
//  ifstream theFile ((fileName).c_str());
//
//  vector<double> row;
//
//  if (theFile.is_open()){
//    while (getline(theFile,line)){ 
//      while (line.length() > line.find(" ")){
//        row.push_back(stod(line.substr(0, line.find(" "))));
//        line = line.substr(line.find(" ") + 1, line.length());
//      }
//      row.push_back(stod(line));      
//      vecToImport.push_back(row[0]);
//      row.erase(row.cbegin(), row.cend());
//    }
//  } else {
//    errorMessage("importdata", " the file "+fileName+" cannot be opened or does not exists.\n");
//  }
//  return vecToImport;
//}
//void importdata1d(string fileName, vector<int> vecToImport){
//
//  string line;
//
//  ifstream theFile ((fileName).c_str());
//
//  vector<double> row;
//
//  if (theFile.is_open()){
//    while (getline(theFile,line)){ 
//      while (line.length() > line.find(" ")){
//        row.push_back(stod(line.substr(0, line.find(" "))));
//        line = line.substr(line.find(" ") + 1, line.length());
//      }
//      row.push_back(stod(line));      
//      vecToImport.push_back(row[0]);
//      row.erase(row.cbegin(), row.cend());
//    }
//  } else {
//    errorMessage("importdata", " the file "+fileName+" cannot be opened or does not exists.\n");
//  }
//}

void deleteRow( int idx, vector<vector<double>>  vec){
  vec[idx] = vec.back();
  vec.pop_back();
}

unsigned int getLine(string filePath, int nLine, string * theLine){
  ifstream theFile(filePath);
  if (theFile.is_open()) {
    for (int i = 0; i < nLine; i++) {
      getline(theFile, *theLine);
    }
  } else {
    return 1;
  }
  theFile.close(); 
  return 0;
}

void exportData(string filePath, vector<vector<int>> dataToExport, string delimiter){
  ofstream File (filePath);     
  for (unsigned int i = 0; i < dataToExport.size(); i++){
    for (unsigned int j = 0; j < dataToExport[0].size(); j++){
      File << dataToExport[i][j] << " ";
    }
    File << endl;
  }
  File.close();
}

void exportData(string filePath, vector<vector<double>> dataToExport, string delimiter, int precision){

  ofstream File (filePath);     

  for (unsigned int i = 0; i < dataToExport.size(); i++){
    for (unsigned int j = 0; j < dataToExport[0].size()-1; j++){
        File << setprecision(precision) << scientific << dataToExport[i][j] << delimiter;
    }
    File << setprecision(precision) << scientific << dataToExport[i][dataToExport[0].size()-1] << endl;
  }
  File.close();
}

void exportData(string filePath, vector<string> dataToExport){

  ofstream File (filePath);     

  for (unsigned int i = 0; i < dataToExport.size(); i++)
      File << dataToExport[i] << endl; 

  File.close();
}

void exportData(string filePath, vector<double> dataToExport){

  ofstream File (filePath);    

  if (File.is_open()){ 
    for (unsigned int i = 0; i < dataToExport.size(); i++)
        File << scientific << dataToExport[i] << endl; 

  } else {
    errorMessage("exportData", "Can't open " + filePath + " for writing.");
    exit(1); 
  }
  File.close();
}

void exportData(vector<double> dataToExport, string filePath){

  ofstream File (filePath);    

  if (File.is_open()){ 
    for (unsigned int i = 0; i < dataToExport.size(); i++)
        File << scientific << dataToExport[i] << endl; 

  } else {
    errorMessage("exportData", "Can't open " + filePath + " for writing.");
    exit(1); 
  }
  File.close();
}

void exportData(string filePath, vector<double> dataToExport, string delimiter, int precision){

  ofstream File (filePath);     

  for (unsigned int i = 0; i < dataToExport.size(); i++){
    File << setprecision(precision) << scientific << dataToExport[i] << delimiter; 
  }

  File.close();
}

void exportData(string filePath, vector<int> dataToExport){

  ofstream File (filePath);     

  for (unsigned int i = 0; i < dataToExport.size(); i++)
      File << to_string( dataToExport[i] ) << endl; 

  File.close();
}

vector<double> vertcat( vector<double> u, vector<double> v ) {
  u.insert( u.end(), v.begin(), v.end() );
  return u;
}

vector<int> findRepeatedIndices(vector<double> u){

  vector<int> redundantElements;
  double eps = 0.00001;

  for (int i = 1; i < u.size() ; i++ ){
    if ( u[i] - u[i-1]  < eps && u[i] - u[i-1] > - eps )
      redundantElements.push_back(i);
  }
  return redundantElements; 
}


void erase(vector<double> & u, vector<int> indices){
  for (int i = 0; i < indices.size(); i++){
    u.erase(u.begin() + indices[i]-i); 
  }
}

void erase(vector<vector<double>> & u, vector<int> indices){
  for (int i = 0; i < indices.size(); i++){
    u.erase(u.begin() + indices[i]-i); 
  }
}

/* erase elements of an STL vector given a vector of increasing indices */
vector<double> eraseElements(vector<double> u, vector<int> indices){
  vector<double> newVector;

  for (int i = 0; i < u.size(); i++)  {
  	if (i == indices[0])
      indices.erase(indices.cbegin());
    else 
      newVector.push_back(u[i]);
  }
  return newVector;
}

/* max and min of a STL vector in Matlab style */
double max(vector<double> u){
  return *(max_element(u.cbegin(), u.cend()));
}
double min(vector<double> u){
  return *(min_element(u.cbegin(), u.cend()));
}

/* find function for stl vectors. Coding:       1: =     2: <       3: >      */
vector<int> findSomething(vector<double> u, double whatToLookFor, int comparison_operator){
  vector<int> indexes;	
  for (int i = 0; i < u.size(); i++) {
  	if ( comparison_operator == 1 ){
      if ( u[i] == whatToLookFor )	
        indexes.push_back(i);       	
    } else if (comparison_operator == 2){
      if ( u[i] < whatToLookFor )	
        indexes.push_back(i);
    } else if (comparison_operator == 3){
      if ( u[i] > whatToLookFor )	
        indexes.push_back(i);
    }
  }
  return indexes;
}

/* find function for stl vectors. Coding:       1: =     2: <       3: >      */
int findSomething(const vector<int> u, int whatToLookFor){
  for (int i = 0; i < u.size(); i++) {
      if ( u[i] == whatToLookFor )	
        return i;
    }
//  errorMessage("findSomething", "Element '" + to_string(whatToLookFor) + "' is not in the list");
  return -1;
}

template<class InputIterator, class T>
  InputIterator find (InputIterator first, InputIterator last, const T& val)
{
  while (first!=last) {
    if (*first==val) return first;
    ++first;
  }
  return last;
}

/* get sub-vector of a stl vector given a vector of incremental indices*/
vector<double> getSubVector(vector<double> u, vector<int> indices){

  vector<double> subVector;
  subVector.resize(indices.size());

  for (int i = 0; i < indices.size(); i++)
    subVector[i] = u[indices[i]];  

  return subVector;
}


/* get some row or column of a STL matrix */
vector<double> getRow(vector<vector<double>> u, int rowId){
  return u[rowId];
}

vector<double> getCol(vector<vector<double>> u, int colId){
  vector<double> theCol;
  theCol.resize(u.size());
  for (int i = 0; i < u.size(); i++){
    theCol[i] = u[i][colId];
  }
  return theCol;
}

/* Returns interpolated value at x from parallel arrays ( xData, yData )
   Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
   boolean argument extrapolate determines behaviour beyond ends of array (if needed)	*/
double interpolate( vector<double> &xData, vector<double> &yData, double x, bool extrapolate)
{
   int size = xData.size();

   int i = 0;                                                                  // find left end of interval for interpolation
   if ( x >= xData[size - 2] ){                                                // special case: beyond right end
      i = size - 2;
   }
   else {
      while ( x > xData[i+1] ) i++;
   }
   double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];      // points on either side (unless beyond ends)
   if ( !extrapolate ){                                                         // if beyond ends of array and not extrapolating
      if ( x < xL ) yR = yL;
      if ( x > xR ) yL = yR;
   }

   double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient

   return yL + dydx * ( x - xL );                                              // linear interpolation
}

/* Matlab like linear interpolation function.
	INPUT:  sampling points stl vectors x y and query points xq
	OUTPUT: interpolated values over query points
*/
vector<double> interp1(vector<double> &x, vector<double> &y, vector<double> &xq){

  vector<double> yq;

  for (auto it = xq.cbegin(); it != xq.cend(); it++){
    yq.push_back( interpolate(x, y, *it) );     
  }

  return yq;
}

/* double vector to string vector */
vector<string> num2Str(vector<double> u){

  vector<string> v;
  v.resize(u.size());

  for (int i = 0; i < u.size(); i++)
    v[i] = to_string(u[i]);

  return v;
}

/* double matrix to string matrix */
vector<vector<string>> num2Str(vector<vector<double>> m){

  vector<vector<string>> strm;
  strm.resize(m.size());

  for (int i = 0; i < m.size(); i++){
  	for (int j = 0 ; j < m[0].size(); j++) {
  	  strm[i].resize(m[0].size());	
      strm[i][j] = to_string(m[i][j]);

    }
  }
  return strm;
}

/* string vector to double vector */
vector<double> str2Num(vector<string> u){

  vector<double> toReturn;
  toReturn.resize(u.size());

  for (int i = 0; i < u.size(); i++)
      toReturn[i] = stod(u[i]);

  return toReturn;

}

/* string matrix to double matrix */
vector<vector<double>> str2Num(vector<vector<string>> m){

  vector<vector<double>> toReturn;
  toReturn.resize(m.size());

  for (int i = 0; i < m.size(); i++){
    for (int j = 0 ; j < m[0].size(); j++){
      toReturn[i].resize(m[0].size());
      toReturn[i][j] = stod(m[i][j]);
    }
  }

  return toReturn;

}

/* horizontal concatenation of vectors matlab-style */
vector<vector<double>> horzcat(vector<double> u, vector<double> v, vector<double> w){
  
  vector<vector<double>> toReturn;
  toReturn.resize(u.size());

  for (int i = 0; i < u.size(); i++) {
  	toReturn[i].push_back(u[i]);
  	toReturn[i].push_back(v[i]);
  	toReturn[i].push_back(w[i]);
  }

  return toReturn;
}

/*
		Methods to handle I/O of binary files for STL vectors. 
*/

/* write binary file from double vector */
void writeBinaryVec(vector<double> u, string fileName){

  ofstream outputFile(fileName,  std::ios::out | std::ios::binary); // saving file

  vector<string> v;

  v = num2Str(u);

  /* save vector size in order to get the read task easier */
  unsigned int vecSize = u.size();
  outputFile.write( (char*)( &vecSize ), sizeof( vecSize ) );

  unsigned int stringLength;
  for (int i = 0; i < u.size(); i++){
    stringLength = v[i].length();
    outputFile.write( (char*)( &stringLength ), sizeof( stringLength ) );
    outputFile.write( v[i].c_str(), v[i].length() );
  }
  outputFile.close();
}

/* write binary file from double matrix */
void writeBinaryMat(vector<vector<double>> m, string fileName){

  ofstream outputFile(fileName,  std::ios::out | std::ios::binary); // saving file

  vector<vector<string>> v;

  v = num2Str(m);

  /* save matrix size in order to get the read task easier */
  unsigned int mRows = m.size();
  unsigned int mCols = m[0].size();
  outputFile.write( (char*)( &mRows ), sizeof( mRows ) );
  outputFile.write( (char*)( &mCols ), sizeof( mCols ) );

  unsigned int elementLength;
  for (int i = 0; i < mRows; i++){
  	for(int j = 0; j < mCols; j++){
      elementLength = v[i][j].length();
      outputFile.write( (char*)( &elementLength ), sizeof( elementLength ) );
      outputFile.write( v[i][j].c_str(), v[i][j].length() );
    }
  }
  outputFile.close();
}

vector<vector<double>> readBinaryMat(string fileName){
  /* read binary file to double STL matrix */
  ifstream inputFile(fileName, std::ios::in | std::ios::binary);	
  /* read vector size */
  unsigned int mRows, mCols;
  inputFile.read( (char*)( &mRows ), sizeof( mRows ) );
  inputFile.read( (char*)( &mCols ), sizeof( mCols ) );
  vector<vector<string>> strVector;
  strVector.resize(mRows);
  unsigned int elementLength;
  for (int i = 0; i < mRows; i++){
  	for (int j = 0; j <  mCols; j++){
      inputFile.read( (char*)( &elementLength ), sizeof( elementLength ) );
      strVector[i].resize(mCols);
      strVector[i][j].resize( elementLength );
      inputFile.read( (char*)strVector[i][j].c_str(), elementLength);
    }
  }
  inputFile.close();
  vector<vector<double>> theVector;
  theVector = str2Num(strVector);
  return theVector;
}


vector<double> readBinaryVec(string fileName){
  /* read binary file to double STL vector */
  ifstream inputFile(fileName, std::ios::in | std::ios::binary);	
  /* read vector size */
  unsigned int vecSize;
  inputFile.read( (char*)( &vecSize ), sizeof( vecSize ) );
  vector<string> strVector;
  strVector.resize(vecSize);
  unsigned int stringLength;
  for (int i = 0; i < vecSize; i++){
    inputFile.read( (char*)( &stringLength ), sizeof( stringLength ) );
    strVector[i].resize( stringLength );
    inputFile.read( (char*)strVector[i].c_str(), stringLength );
  }
  inputFile.close();
  vector<double> theVector;
  theVector = str2Num(strVector);
  return theVector;
}

double determinant44(const vector<vector<double>> & m){
  vector<vector<double>> m1(3);
  m1[0] = m[1];
  m1[0].erase(m1[0].begin()+0);
  m1[1] = m[2];
  m1[1].erase(m1[1].begin()+0);
  m1[2] = m[3];
  m1[2].erase(m1[2].begin()+0);

  vector<vector<double>> m2(3);
  m2[0] = m[1];
  m2[0].erase(m2[0].begin()+1);
  m2[1] = m[2];
  m2[1].erase(m2[1].begin()+1);
  m2[2] = m[3];
  m2[2].erase(m2[2].begin()+1);

  vector<vector<double>> m3(3);
  m3[0] = m[1];
  m3[0].erase(m3[0].begin()+2);
  m3[1] = m[2];
  m3[1].erase(m3[1].begin()+2);
  m3[2] = m[3];
  m3[2].erase(m3[2].begin()+2);

  vector<vector<double>> m4(3);
  m4[0] = m[1];
  m4[0].erase(m4[0].begin()+3);
  m4[1] = m[2];
  m4[1].erase(m4[1].begin()+3);
  m4[2] = m[3];
  m4[2].erase(m4[2].begin()+3);

  return m[0][0] * determinant33(m1) - m[0][1] * determinant33(m2) 
       + m[0][2] * determinant33(m3) - m[0][3] * determinant33(m4);
}

double determinant33(const vector<vector<double>> & m){
  return m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2]) - m[0][1]*(m[1][0]*m[2][2] - m[2][0]*m[1][2]) + m[0][2]*(m[1][0]*m[2][1] - m[2][0]*m[1][1]);
}

double determinant22(const vector<vector<double>> & m){
  return m[0][0] * m[1][1] - m[1][0]*m[0][1];
}

double determinant(const vector<vector<double>> & m){
  assert(m.size() == m[0].size()); /* only for 3x3 or 2x2 matrices */
  if (m.size() == 1){
    return m[0][0]; 
  } else if (m.size() == 2){
    return determinant22(m);
  } else if (m.size() == 3){
    return determinant33(m);
  } else if (m.size() == 4) {
    return determinant44(m);
  } else {
    cout << "Trying to compute determinant of a non supported size." << endl;
    exit(1);
  }
}

vector<vector<double>> invert(const vector<vector<double>> & m){
  assert(m.size() == 3 || m.size() == 2 || m.size() == 1); 
  if (m.size() == 3){
    assert(m[0].size() == 3); /* only for 3x3 matrices */
    assert(m[1].size() == 3); /* only for 3x3 matrices */
    assert(m[2].size() == 3); /* only for 3x3 matrices */
    /* adjugate of m */
    vector<vector<double>> adj_m; 
    adj_m.resize(3);
    adj_m[0].resize(3);
    adj_m[0][0] =    m[1][1]*m[2][2] - m[2][1]*m[1][2];
    adj_m[0][1] = - (m[1][0]*m[2][2] - m[2][0]*m[1][2]);
    adj_m[0][2] =    m[1][0]*m[2][1] - m[2][0]*m[1][1];
    adj_m[1].resize(3);
    adj_m[1][0] = - (m[0][1]*m[2][2] - m[2][1]*m[0][2]);
    adj_m[1][1] =    m[0][0]*m[2][2] - m[2][0]*m[0][2];
    adj_m[1][2] = - (m[0][0]*m[2][1] - m[2][0]*m[0][1]);
    adj_m[2].resize(3);
    adj_m[2][0] =    m[0][1]*m[1][2] - m[1][1]*m[0][2];
    adj_m[2][1] = - (m[0][0]*m[1][2] - m[1][0]*m[0][2]);
    adj_m[2][2] =    m[0][0]*m[1][1] - m[1][0]*m[0][1];
    return 1.0/determinant(m)*transpose(adj_m);

  } else if (m.size() == 2){
    assert(m[0].size() == 2); /* only for 2x2 matrices */
    assert(m[1].size() == 2); /* only for 2x2 matrices */
    /* Adjugate of m */
    vector<vector<double>> adj_m; 
    adj_m.resize(2);
    adj_m[0].resize(2);
    adj_m[0][0] =   m[1][1];
    adj_m[0][1] = - m[0][1];
    adj_m[1].resize(2);
    adj_m[1][1] =   m[0][0];
    adj_m[1][0] = - m[1][0];

    return 1.0 / determinant(m) * adj_m; 
  } else {
    vector<vector<double>> inv(1);
    inv[0].push_back(1.0/m[0][0]); 
    return inv;
  }
}

vector<vector<double>> transpose(const vector<vector<double>> & m){
  vector<vector<double>> m_t(m[0].size());
  for(int i = 0; i < m_t.size(); i++){
    m_t[i].resize(m.size());
    for(int j = 0; j < m_t[i].size(); j++){
      m_t[i][j] = m[j][i];
    }
  }
  return m_t;
}

vector<vector<double>> transpose(const vector<double> & m){
  vector<vector<double>> m_t(1);
  m_t[0].resize(m.size());
  for(int j = 0; j < m.size(); j++){
    m_t[0][j] = m[j];
  }
  return m_t;
}

vector<vector<int>> transpose(const vector<vector<int>> & m){
  vector<vector<int>> m_t(m[0].size());
  for(int i = 0; i < m.size(); i++){
    m_t[i].resize(m.size());
    for(int j = 0; j < m[i].size(); j++){
      m_t[i][j] = m[j][i];
    }
  }
  return m_t;
}

double trace(const vector<vector<double>> & m){
  assert(m.size() == m[0].size());
  double trace = 0.0;
  for (int i = 0; i < m.size(); i++){
    trace += m[i][i];
  }
  return trace;  
}

vector<double> arcsin(vector<double> u){
  vector<double> v(u.size());
  for (int i=0; i<u.size(); i++){
    v[i] = asin(u[i]);
  }
  return v;
}

vector<double> arccos(vector<double> u){
  vector<double> v(u.size());
  for (int i=0; i<u.size(); i++){
    v[i] = asin(u[i]);
  }
  return v;
}

double average(vector<double> const& u){
  if(u.empty()){
    return 0;
  }
  double count = static_cast<double>(u.size());
  return sum(u) / count;
}
