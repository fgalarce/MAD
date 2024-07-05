#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <fstream>
#include "tools_exec.h"
#include "parameter_recogn.h"
#include "mini/ini.h"
#include "mm_file.hpp"


std::tuple<mINI::INIStructure, std::map<std::string, bool>> f_LoadPar(std::string par_name) {
	std::map<std::string, bool> aBool;
	mINI::INIStructure current_par_ini;
	std::vector<std::string> secc;
	std::vector<std::string> keys;
	std::vector<std::string> pars;
	std::vector<int> param_numb;
	
	//Load the par into an ini structure
	mINI::INIFile file(par_name);
	file.read(current_par_ini);

	//Counts the number of sections
	//int secc_numb = current_par_ini.size();

	//Get the par info
	std::tuple< std::vector<std::string>,std::vector<std::string>,std::vector<std::string>,std::vector<int> > par_info=GetParInfo(current_par_ini);
	keys = std::get<0>(par_info);
	pars = std::get<1>(par_info);
	secc = std::get<2>(par_info);
	param_numb = std::get<3>(par_info);
	
	//Sets aBool
	std::map<std::string, int> ParamManifest=GetParamManifest();
	//std::map<std::string, bool> aBool=GetDefaultBools(keys, pars, ParamManifest);
	aBool=GetDefaultBools(keys, pars, ParamManifest);
	//std::cout << aBool["monitorKSP"];
	//std::cout << aBool["monitorksp"];
	
	std::tuple<mINI::INIStructure, std::map<std::string, bool>> LoadedPar;
	LoadedPar = std::make_tuple(current_par_ini,aBool);
	return LoadedPar;
}


mINI::INIStructure CreatePar(std::vector<std::string> secc, std::vector<std::string> keys, std::vector<std::string> pars, std::vector<int> param_numb) {
	//Creates an INIStructure based on the data given
	mINI::INIStructure new_par;
	int current_par = 0;
	for (int i=0;i<secc.size();i++) {
		for (int j=0; j<param_numb[i]; j++) {
			new_par[secc[i]][keys[current_par]] = pars[current_par];
			current_par++;
		}
	}
	return new_par;
}


void f_SavePar(std::string par_name, mINI::INIStructure &current_par_ini) {
	//Saves new data.
	std::tuple< std::vector<std::string>,std::vector<std::string>,std::vector<std::string>,std::vector<int> > par_info;
	par_info = GetParInfo(current_par_ini);
	std::vector<std::string> secc = std::get<2>(par_info);
	std::vector<std::string> keys = std::get<0>(par_info);
	std::vector<std::string> pars = std::get<1>(par_info);
	std::vector<int> param_numb = std::get<3>(par_info);

	std::cout << "---Saving--- \n";
	int k = 0;
	for (int i=0; i<secc.size(); i++) {
		std::cout << secc[i] << "\n";
		for (int j=0; j<param_numb[i]; j++) {
			std::cout << k << keys[k] << " " << pars[k] << "\n";
			current_par_ini[secc[i]][keys[k]] = pars[k];
			k++;
		}
	}
	mINI::INIFile file(par_name);
	file.generate(current_par_ini);
	std::cout << "---Saved--- \n";
}


void f_RunMad(int np, std::string exe_name, mINI::INIStructure current_par_ini) {
	//Create a temporal ".par" based on the changes to the loaded par.
	std::string par_name = ".par";
	f_SavePar(par_name, current_par_ini);
	
	//Runs MAD
	std::string run_mad = "$PETSC_ROOT/build/bin/mpirun -np " + std::to_string(np) + " ./" + exe_name + " " + par_name + " &";
	std::system(run_mad.c_str());
}

void f_RunMadness() {
	//This is not illegal, right?
	//This will compile***
}

void f_RunParaview(std::string paraview_exec, std::string results_folder, std::string patient_name) {
	std::system("echo Opening Paraview...");
	std::string open_parav=paraview_exec + " " + results_folder + patient_name +".case";
	std::system(("echo "+open_parav).c_str());
	std::system(open_parav.c_str());
}


double f_GetDoubleFromBin(std::string filename, int id) {
	// instruct the kernel that we will read the content
    // of the file sequentially
    int advice = mm::advice::sequential;
    // read the stream as doubles
    mm::file_source<double> fin(filename, advice);
    auto const* data = fin.data();
	double value = data[id];
    fin.close();
	return value;
}


void f_WriteDoubleToBin(std::string file, int id, double value, bool bClose) {
	//static bool bFirstTime = true;
	//if (bFirstTime) {
		//mm::file_sink<double> fout(file);
		//bFirstTime = false;
	//}
	//"feedback.bin" is created in gui/src/main.cpp.
	std::ifstream file_iq;
	file_iq.open("feedback.bin");
	
	if (file_iq) {
		mm::file_sink<double> fout(file);
		auto* data = fout.data();
		data[id] = value;
		//for (int i=0; i<fout.size(); ++i) {
		//	data[i] = value;
		//}
		
		if (bClose) {
			fout.close();
		}
	}
}

//sadaddss