#ifndef TOOLS_MISC
#define TOOLS_MISC

#include <string>
#include <vector>

//Get the index of a string in a vector string
int getIndex(std::vector<std::string> v, std::string K);
//String to Bool
bool stb (std::string sbool);
//Bool to String
std::string bts (bool bbool);
//Returns the name of the executable
std::string GetExeName();

#endif