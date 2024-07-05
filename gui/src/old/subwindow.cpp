//#include "subwindow.hpp"
#include <cstring>
#include <iostream>
#include <vector>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "mini/ini.h"

//create a structure ("ini") that will hold the "par" data
mINI::INIStructure ini;


void load_ini(char par_name[128]) {
    //std::string testt = "par";
    mINI::INIFile file(par_name);
    file.read(ini);
    std::cout << par_name;
}

bool show_welcome = true;
bool show_par = false;
bool btest = true;
char par_name[128] = "par";
char par_name1[128] = "par1";
std::vector<char*> teststr(8,"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa");
std::vector<char*> pars(128,"qwertyuiop");


void w_welcome() {
    //Displays a "welcome" text
    ImGui::Text("Welcome to Multi-physics simulAtions for engineering and Data assimilation\n(or MAD for short).\nTry loading a project by choosing its parameters file.");
}


void w_par() {
    //int i = 0;
    //HAY VARIABLES QUE NO SE USAN EN EL CICLO
    /*
    for (auto const& it : ini) {
	    auto const& section = it.first;
	    auto const& collection = it.second;
        ImGui::Text(("---"+section+"---").c_str());
	    for (auto const& it2 : collection) {
		    auto const& key = it2.first;
		    auto const& value = it2.second;
            ImGui::Text((key+" ").c_str());
            ImGui::SameLine();
            ImGui::InputText(" ", pars[i], IM_ARRAYSIZE(pars[i]));
            i++;
	    }
    }
    */
    ImGui::InputText(" ", pars[0], IM_ARRAYSIZE(pars[0]));
    ImGui::InputText(" ", pars[1], IM_ARRAYSIZE(pars[1]));

    if (ImGui::Button("Run simulation")) {
        std::string par_string_name = par_name;
        std::string run_mad="/home/septem/Documentos/MAD/petsc/build/bin/mpirun -np 4 ./ns.exe "+par_string_name;
        std::system(run_mad.c_str());
    }
    if (ImGui::Button("See results in PARAVIEW")) {
        std::system("echo Opening Paraview...");
        //This should be setted up in a configuration file...
        std::string open_parav="/home/septem/paraview/bin/paraview "+ ini["IO"]["dirResults"]+ini["IO"]["patientName"]+".case";
        std::system(("echo "+open_parav).c_str());
        std::system(open_parav.c_str()); //+ ini["IO"]["dirResults"]
    }
}


void w_default() {
    //Default window elements. Also controls what is shown at real time.
    ImGui::Begin("Simulation parameters");
    //ImGui::InputText("<-", par_name, IM_ARRAYSIZE(par_name));
    //par_name1 = ImGui::InputText("<-", par_name1, IM_ARRAYSIZE(par_name1));
    static std::vector<char> buf(2,"aaaaaaaa");
    ImGui::InputText("default1",     buf[0], 32);
    ImGui::InputText("default2",     buf[1], 32);
    ImGui::SameLine();
    if (ImGui::Button("Load par")) {
        show_par = true;
        //DEBERÌA PRIMERO COMPROBAR SI ES VÀLIDO
        //Load the ini
        load_ini(par_name);
        //Completes the "pars" array
        int par_numb = 0;
        for (auto const& it : ini) {
	        auto const& section = it.first;
	        auto const& collection = it.second;
	        for (auto const& it2 : collection) {
    		    auto const& key = it2.first;
		        auto const& value = it2.second;
                //ONLY ACEPTS 128 PARAMETERS
                char *cstr = new char[128];
                strcpy(cstr, value.c_str());
                pars[par_numb] = cstr;
                par_numb++;
	        }
        }
    }
    if (show_welcome) {
        w_welcome();
    }
    if (show_par) {
        //Display the contents of the ini file
        w_par();
    }

    ImGui::End();
}

/*
        // We use a Begin/End pair to create a named window.
        //{
            static float f = 0.0f;
            static int counter = 0;
            static char str0[128] = "Hello, world!";
            //static float timeStep = par.timeStep();
            //static float nbIterations =par.nbIterations();
            //static float amplitude = par.amplitude();
            //static int outlet = par.outlet();

            ImGui::Begin("Simulation parameters");

            if (ImGui::Button("Save par"))
                counter++;
            ImGui::SameLine();
            if (ImGui::Button("Load par"))
                counter++;
            ImGui::Text("------------");
            test1();
            ImGui::Checkbox("3D", &td_sim);

            for (auto const& it : ini) {
	            auto const& section = it.first;
	            auto const& collection = it.second;
                ImGui::Text(("---"+section+"---").c_str());
	            for (auto const& it2 : collection) {
		            auto const& key = it2.first;
		            auto const& value = it2.second;
                    //char str0[128] = value; //HERE
                    ImGui::Text((key+" "+value).c_str());
                    ImGui::SameLine();
                    ImGui::InputText("<-", str0, IM_ARRAYSIZE(str0));
	            }
            }


            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
            ImGui::End();
        //}
*/