#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <fstream>


/* General use functions */


int getIndex(std::vector<std::string> v, std::string K) {
	//Return the index of certain string inside a vector. -1 means not found.
    auto it = find(v.begin(), v.end(), K); 

    if (it != v.end())  
    { 
        int index = it - v.begin(); 
        return index; 
    } 
    else { 
        return -1; 
    } 
}


bool stb (std::string sbool) {
	//Input a string, returns a bool
	if (sbool=="false" || sbool=="0" || sbool=="no" || sbool=="False") {
		return false;
	} else {
		return true;
	  }
}


std::string bts (bool bbool) {
	if (bbool) {
		return "true";
	} else {
		return "false";
	}
}

std::string GetExeName() {
	std::ifstream file("makefile");
	std::string line;
	std::string exe_name;
	if(file.is_open()) {
		while(std::getline(file, line)) {
			if(line.substr(0,12)=="project_name") {
				exe_name= line.substr(13)+".exe";
				break;
			}
		}
	}
	return exe_name;
}