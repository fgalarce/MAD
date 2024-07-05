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

#include<tools.hpp>

string wildcard(string i, string _i){
  return wildcard(stoi(i));
}
string wildcard(size_t i, string _i){
  string iteration;
  if (i < 10)
    iteration = _i + _i + _i + _i + to_string(i);
  else if (i < 100)
    iteration = _i + _i + _i + to_string(i);
  else if (i < 1000)
    iteration = _i + _i + to_string(i);
  else if (i < 10000)
    iteration = _i + to_string(i);
  else if (i < 100000)
    iteration = to_string(i);
  return iteration;
}

string wildcard(int i, string _i){
  string iteration;
  if (i < 10)
    iteration = _i + _i + _i + _i + to_string(i);
  else if (i < 100)
    iteration = _i + _i + _i + to_string(i);
  else if (i < 1000)
    iteration = _i + _i + to_string(i);
  else if (i < 10000)
    iteration = _i + to_string(i);
  else if (i < 100000)
    iteration = to_string(i);
  return iteration;
}

string wildcard8(size_t i, string _i){
  string iteration;
  if (i < 10)
    iteration = _i + _i + _i + _i + _i + _i + _i + to_string(i);
  else if (i < 100)
    iteration = _i + _i + _i + _i + _i + _i + to_string(i);
  else if (i < 1000)
    iteration = _i + _i + _i + _i + _i + to_string(i);
  else if (i < 10000)
    iteration = _i + _i + _i + _i + to_string(i);
  else if (i < 100000)
    iteration = _i + _i + _i + to_string(i);
  else if (i < 1000000)
    iteration = _i + _i +  to_string(i);
  else if (i < 10000000)
    iteration = _i + to_string(i);
  else if (i < 100000000)
    iteration = to_string(i);
  return iteration;
}

string wildcard8(int i, string _i){
  string iteration;
  if (i < 10)
    iteration = _i + _i + _i + _i + _i + _i + _i + to_string(i);
  else if (i < 100)
    iteration = _i + _i + _i + _i + _i + _i + to_string(i);
  else if (i < 1000)
    iteration = _i + _i + _i + _i + _i + to_string(i);
  else if (i < 10000)
    iteration = _i + _i + _i + _i + to_string(i);
  else if (i < 100000)
    iteration = _i + _i + _i + to_string(i);
  else if (i < 1000000)
    iteration = _i + _i +  to_string(i);
  else if (i < 10000000)
    iteration = _i + to_string(i);
  else if (i < 100000000)
    iteration = to_string(i);
  return iteration;
}

string get_filename_extension(string filename) {

  string filename_inverted = filename;
  
  reverse(filename_inverted.begin(), filename_inverted.end());

  string extension = filename_inverted.substr(0, filename_inverted.find("."));

  reverse(extension.begin(), extension.end());
  return extension;

}
//
///* search for line top and replace it by bottom */
//void fileReplace(string top, string bottom, string fileName){
//  ifstream
//
//}
//
//def fileReplace(top, bottom, fileName):
//  tempFile =  open(fileName, 'r+' )
//  for line in fileinput.input( fileName ):
//      tempFile.write( line.replace(top, bottom)  )
//  tempFile.close()
   
/*!
 Splits the string str by any char appears in delim, not including this char.
 For instance, the call SplitString( "abcdefghabcd", "ade", v);
 will put in v the strings "bc" "fgh" "bc".
 */
void splitString(string str, vector<string>& vect, string delimiter) {
  vect.clear();

  string tmp;
  string::size_type index_beg, index_end;

  index_beg = str.find_first_not_of(delimiter);

  while (index_beg != string::npos) {
    index_end = str.find_first_of(delimiter, index_beg);
    convert(str.substr(index_beg, index_end == string::npos ?
                       string::npos : (index_end - index_beg)), tmp);
    vect.push_back(tmp);
    index_beg = str.find_first_not_of(delimiter, index_end);
  }
}

void convert(const std::string& s, std::string& out) {
  out = s;
}

string getParameter(string parameter, string filename){
  
  check_existence(filename);

  string data;
  bool found_data = false;

  ifstream dataFile(filename);
  if (dataFile.is_open()){
    string line;
    while(getline(dataFile,line)){
      if (line.substr(0, line.find("=")) == parameter){
        found_data = true;
        data = line.substr(line.find("=")+1, line.back());
        return data;
      }
    }
  }

  dataFile.close();

  if (!found_data)
    MDMA_error("Imposible to find variable " + parameter + " in: " + filename);
  return "ERROR";
}

void cat(string file){
  ifstream parametersFile(file);
  string line; 
  cout << endl;
  while (getline(parametersFile,line))
    cout << line << endl;
}

void mkdir(string folder){
  if( system(("mkdir -p " + folder).c_str()) < 0) {
    errorMessage("mkdir", "problems creating folder " + folder);
  }
}

void check_existence(string file){
  ifstream exist_test(file);
  if(!exist_test.good()){
    errorMessage("check_existence", "problem reading " + file);
  }
}

void errorMessage(string functionName, string theError){
  int world_rank;  
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank == 0) cout << "\n - - - - - - - - - "<< "\033[1;31mMAD ERROR\033[0m" << " - - - - - - - - - " << endl;
  if (world_rank == 0) cout << "In function: " << functionName << endl;
  if (world_rank == 0) cout << "Cause: " << theError;
  if (world_rank == 0) cout << "\n - - - - - - - - - "<< "\033[1;31mMAD ERROR\033[0m" << " - - - - - - - - - " << endl;
  exit(1);
}

void print(string message){
  int world_rank;  
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank == 0) cout << "\n - - - - - - - - - "<< "\033[1;31mMAD\033[0m" << " - - - - - - - - - " << endl;
  if (world_rank == 0) cout << message;
  if (world_rank == 0) cout << "\n - - - - - - - - - "<< "\033[1;31mMAD\033[0m" << " - - - - - - - - - " << endl;
}

void MDMA_error(string theError){
  cout << "\n - - - - - - - - - "<< "\033[1;31mERROR\033[0m" << " - - - - - - - - - " << endl;
  cout << "\t" << theError;
  cout << "\n - - - - - - - - - "<< "\033[1;31mERROR\033[0m" << " - - - - - - - - - " << endl;
  exit(1);
}

string get_env_var(string const & key ) {                                 
  char * val;                                                                        
  val = getenv( key.c_str() );                                                       
  std::string retval = "";                                                           
  if (val != NULL) {                                                                 
    retval = val;                                                                    
  }                                                                                  
  return retval;                                                                        
}   

void loop_cout(int iter, int n, string message){
  int world_rank;  
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (iter % (n / 5) == 0){
    if (world_rank == 0) cout << message << endl;
  }
}

void loading_bar(int iter, int n){
  int world_rank;  
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (iter % (n / 5) == 0){
    if (world_rank == 0) cout << "- ";
  }
}
