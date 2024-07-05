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

#ifndef TOOLS
#define TOOLS 

#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <string>
#include <chrono>
#include <ctime>
#include <fstream>
#include <slepc.h>
#include <stdlib.h>

using namespace std;

class Timer {
private:
    // Use the high_resolution_clock for the best precision
    std::chrono::high_resolution_clock::time_point start_time;

public:
    // Constructor
    Timer() {
        // Automatically start the timer upon creation
        tic();
    }

    // Start or restart the timer
    void tic() {
        start_time = std::chrono::high_resolution_clock::now();
    }

    // Stop the timer and print the elapsed time in seconds
    double toc() {
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
        double seconds = duration.count() / 1000000.0;
        return seconds;
    }

};



/* get/set parameter in data files */
string getParameter(string parameter, string filename);

string wildcard(string i, string _i = "0");
string wildcard(size_t i, string _i = "0");
string wildcard8(size_t i, string _i = " ");
string wildcard(int i, string _i = "0");
string wildcard8(int i, string _i = " ");

string get_filename_extension(string filename);

/* function taken from Felisce */
void splitString(string str, vector<string>& vect, string delimiter = " ");

/* function taken from Felisce */
void convert(const string& s, string& out);

string getFromDataFile(string dataFilePath, string variableName);

/* bash like functions */
void mkdir(string folder);
void cat(string file);

void check_existence(string file);

void print(string message);
void errorMessage(string functionName, string error);
void MDMA_error(string theError);
void loop_cout(int iter, int n, string message);
void loading_bar(int iter, int n);
string get_env_var(string const & key);
#endif
