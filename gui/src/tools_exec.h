#ifndef TOOLS_EXEC
#define TOOLS_EXEC

#include <string>
#include <vector>
#include <map>
#include "mini/ini.h"
#include "mm_file.hpp"


//Load a par file into a mINI::INIStructure variable
std::tuple<mINI::INIStructure, std::map<std::string, bool>> f_LoadPar(std::string par_name);
//Saves a par file
void f_SavePar(std::string par_name, mINI::INIStructure &current_par_ini); 
mINI::INIStructure CreatePar(std::vector<std::string> secc, std::vector<std::string> keys, std::vector<std::string> pars, std::vector<int> param_numb);
//Run a simulation in MAD
void f_RunMad(int np, std::string exe_name, mINI::INIStructure current_par_ini);
//Display the results in PARAVIEW
void f_RunParaview(std::string paraview_exec, std::string results_folder, std::string patient_name);
//Returns a double from a binary file
double f_GetDoubleFromBin(std::string filename, int id);
//Writes a double to a binary file
void f_WriteDoubleToBin(std::string file, int id, double value, bool bClose);

#endif