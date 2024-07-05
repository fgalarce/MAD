#ifndef PARAMETER_RECOGN
#define PARAMETER_RECOGN


#include <string>
#include <vector>
#include <map>
#include "mini/ini.h"

//Returns the ParamManifest (defined in parameter_recogn.cpp) as a map.
std::map<std::string, int> GetParamManifest();

std::map<std::string, std::vector<std::string>> GetaDropdown();

//Returns the type (according to the ParamType enum)
int GetParamType (std::string ParamKey, std::map<std::string, int> ParamManifest);

std::map<std::string, bool> GetDefaultBools (std::vector<std::string> keys, std::vector<std::string> pars, std::map<std::string, int> ParamManifest);

//Returns keys, values and sections y a tuple of string vectors
std::tuple< std::vector<std::string>,std::vector<std::string>,std::vector<std::string>,std::vector<int> > GetParInfo(mINI::INIStructure ini);

//Returns a manifest of aliases
std::tuple< std::map<std::string, std::string>, std::map<std::string, std::string> > GetParamAlias();


#endif