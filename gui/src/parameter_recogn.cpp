#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <tuple>
#include "mini/ini.h"
#include "parameter_recogn.h"
#include "tools_misc.h"


enum ParamType{
	tDef=0, //0 (Default text field)
	tBool, //1 (Checkbox)
	tInt, //2
	tFloat, //3
	tDropDown=4 //4 (Dropdown menu)
	};


std::map<std::string, int> GetParamManifest() {
	std::map<std::string, int> ParamManifest;
	//By now, only works 0, 1 and 4.
	//ParamManifest[""] = "";
	ParamManifest["hemodynamics"] = tBool;
	ParamManifest["backflowStab"] = tBool;
	ParamManifest["writeNonLinearIterations"] = tBool;
	ParamManifest["solver"] = tDropDown;
	ParamManifest["useModifiedGramSchmidt"] = tBool;
	ParamManifest["reuse_preconditioner"] = tBool;
	ParamManifest["use_solution_as_guess_KSP"] = tBool;
	ParamManifest["monitorKSP"] = tBool;
	ParamManifest["preconditioner"] = tDropDown;
	
	return ParamManifest;
}


std::tuple< std::map<std::string, std::string>, std::map<std::string, std::string> > GetParamAlias() {
	//Returns the alias and the tooltip that will appear in the UI for certain parameter.
	//If not found, the raw name (as in the par file) will appear
	std::map<std::string, std::string> ParamAlias;
	std::map<std::string, std::string> ToolTips;
	ParamAlias["backflowStab"] = "Stabilize backflow?";
	ToolTips["backflowStab"] = "Defines... PLACEHOLDER TEXT";
	ParamAlias["reuse_preconditioner"] = "Reuse preconditioner?";
	ToolTips["reuse_preconditioner"] = "States that...";
	std::tuple< std::map<std::string, std::string>, std::map<std::string, std::string> > Aliases = std::make_tuple(ParamAlias, ToolTips);
	return Aliases;
}


std::map<std::string, std::vector<std::string>> GetaDropdown() {
	std::map<std::string, std::vector<std::string>> aDropdown;
	//aDropdown["xx"] stores the variable names that MAD can recognize
	//aDropdown["xx_alias"] stores the names that the UI shows
	//aDropdown[""] = {"", "", ...};
	aDropdown["solver_alias"] = {"Conjugate Gradient","Generalized MinRES", "preonly"};
	aDropdown["solver"] = {"cg", "gmres", "preonly"};
	aDropdown["preconditioner_alias"] = {"Additive Schwarz", "Multi-grid"};
	aDropdown["preconditioner"] = {"asm", "mg"};
	
	return aDropdown;
}


std::map<std::string, bool> GetDefaultBools (std::vector<std::string> keys, std::vector<std::string> pars, std::map<std::string, int> ParamManifest) {
	std::map<std::string, bool> aBool;
	//std::cout << "GetDefaultBools()"<< "\n";
	for (int i=0; i<keys.size(); i++) {
		if (ParamManifest.count(keys[i]) && ParamManifest[keys[i]]==1) {
			//std::cout << "GetDefaultBools()"<< "\n";
			aBool[keys[i]] = stb(pars[i]);
			//std::cout << pars[i] << "\n";
			//std::cout << aBool[keys[i]] << "\n";
		}
	}
	
	return aBool;
}


int GetParamType (std::string ParamKey, std::map<std::string, int> ParamManifest) {
	//Returns the key's type as an index of the ParamType enum
	if(ParamManifest.count(ParamKey)) {
		//If the type is known, this will happen
		int ParamType = ParamManifest[ParamKey];
		return ParamType;
	} else { 
		return 0;
	  }
}


std::tuple< std::vector<std::string>,std::vector<std::string>,std::vector<std::string>,std::vector<int> > GetParInfo(mINI::INIStructure ini) {
	//Returns a tuple with the keys, values and sections as string vectors
	std::vector<std::string> pars; //Vector with the parameters
	std::vector<std::string> keys; //Vector with the parameter keys
	std::vector<std::string> secc; //Vector with the sections
	std::vector<int> param_numb; //Vector with the number of parameters per section
	
	for (auto const& it : ini) {
		auto const& section = it.first;
		secc.push_back(section);
		auto const& collection = it.second;
		for (auto const& it2 : collection) {
			auto const& key = it2.first;
			auto const& value = it2.second;
			keys.push_back(key);
			pars.push_back(value);
		}
	}
	for (int i=0; i<ini.size(); i++) {
		int num = ini[secc[i]].size();
		param_numb.push_back(num);
	}
	
	std::tuple< std::vector<std::string>,std::vector<std::string>,std::vector<std::string>,std::vector<int> > inituple(keys,pars,secc,param_numb);
	
	return inituple;
}
