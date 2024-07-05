#include <iostream>
#include <thread>
#include <string>
#include <map>
#include <vector>
#include <tuple>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
//#include "misc/cpp/imgui_stdlib.h"
#include "misc/cpp/imgui_stdlib.cpp"
#include "mini/ini.h"
#include "parameter_recogn.h"
#include "subwindow.h"
#include "tools_exec.h"
#include "tools_misc.h"
#include "implot.h"
#include <math.h>

enum window_types{
	par_editor=0,
	graph_test,
	plot_demo
};

//The par file that is loaded and can be saved
static mINI::INIStructure current_par_ini;
//Temporal par file that will be launched
static mINI::INIStructure temp_par_ini;
//keys, pars y secc
static std::vector<std::string> secc;
static std::vector<std::string> keys;
static std::vector<std::string> pars;
static std::vector<int> param_numb;
static std::vector<int> item_current_idx(128,0);


void show_subwindow(int win_type) {
	switch (win_type) {
		case par_editor:
		{
			display_ParEditor();
			break;
		}
		case graph_test:
		{
			display_GraphTest();
			break;
		}
		case plot_demo:
		{
			display_DemoPlot();
			break;
		}
	}
}


static std::map<std::string,int> ParamManifest = GetParamManifest();
static std::tuple< std::map<std::string, std::string>, std::map<std::string, std::string> > Aliases = GetParamAlias();
static std::map<std::string, std::string> ParamAlias = std::get<0>(Aliases);
static std::map<std::string, std::string> ToolTips = std::get<1>(Aliases);
static std::map<std::string,std::vector<std::string>> aDropdown = GetaDropdown();
static std::map<std::string, bool> aBool;


static void HelpMarker(const char* desc) {
    ImGui::TextDisabled("(?)");
    if (ImGui::BeginItemTooltip())
    {
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
        ImGui::TextUnformatted(desc);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}


void DisplayParam(std::string ParamKey, std::string &ParamRealValue,std::string IdInput, int ParamNumber) {
	int ParamType = GetParamType(ParamKey, ParamManifest);
	//std::cout << ParamKey;
	if (ToolTips.count(ParamKey)) {
		//Create hoverable <?> with a tooltip (if defined)
		std::string tooltip = ToolTips[ParamKey];
		HelpMarker(tooltip.c_str());
		ImGui::SameLine();
	}
	switch(ParamType) {
		case 1: //tBool
		{
			ImGui::Checkbox(IdInput.c_str(), &aBool[ParamKey]);
			pars[ParamNumber] = bts(aBool[ParamKey]);
			break;
		}
		case 20: //tInt
		{
			ImGui::Text("PLACEHOLDER");
			std::string int_id="##int"+ParamKey;
			//ImGui::InputInt(int_id.c_str(), &x)
			break;
		}
		case 30: //tFloat
		{
			ImGui::Text("PLACEHOLDER");
			break;
		}
		case 4: //tDropDown
		{
			std::string combo_id = "##combo"+ParamKey;
			//item_current_idx debe ser diferente para cada dropdown
			//static int item_current_idx = getIndex(aDropdown["solver"], ParamRealValue);
			static std::vector<int> item_current_idx(128,0);
			item_current_idx[ParamNumber] = getIndex(aDropdown[ParamKey], ParamRealValue);
			if (item_current_idx[ParamNumber]<0) {
				item_current_idx[ParamNumber]=0;
			}
			//std::cout << "\n" << ParamRealValue << ParamNumber <<item_current_idx[ParamNumber];
			std::string combo_preview_value = aDropdown[ParamKey+"_alias"][item_current_idx[ParamNumber]];
			
			if (ImGui::BeginCombo(combo_id.c_str(), combo_preview_value.c_str())) {
				for (int n = 0; n < aDropdown[ParamKey].size(); n++) {
					const bool is_selected = (item_current_idx[ParamNumber] == n);
					if (ImGui::Selectable(aDropdown[ParamKey+"_alias"][n].c_str(), is_selected))
						item_current_idx[ParamNumber] = n;

					// Set the initial focus when opening the combo (scrolling + keyboard navigation focus)
					if (is_selected)
						ImGui::SetItemDefaultFocus();
				}
				ImGui::EndCombo();
			}
			pars[ParamNumber] = aDropdown[ParamKey][item_current_idx[ParamNumber]];
			break;
		}
		default:
			ImGui::InputText(IdInput.c_str(), &pars[ParamNumber]);
			break;
	}
}


void display_ParEditor() {
	static char par_name[32]="par"; //Imported filename
	static char par_name_to_export[32]="par"; //Filename to export
	static bool show_pars = false;  //Shows the parameters only if a file is loaded
	static std::vector<int> param_numb; //(32,0);  //Number of parameters in every section
	static bool bTryToLoad = true;
	
	ImGui::Begin("PAR editor");
	ImGui::Text("Welcome to Multi-physics simulAtions for engineering and Data assimilation\n(or MAD for short).\nTry loading a project by choosing its parameters file.");
	ImGui::InputText("##Input", par_name, 32);
	
	if (bTryToLoad) {
		//Loads the par into a INIStructure
		std::tuple<mINI::INIStructure, std::map<std::string, bool>> LoadedPar;
		LoadedPar = f_LoadPar(par_name);
		current_par_ini = std::get<0>(LoadedPar);
		aBool = std::get<1>(LoadedPar);
		//Sets keys, pars and secc
		std::tuple< std::vector<std::string>,std::vector<std::string>,std::vector<std::string>,std::vector<int> > par_info;
		par_info = GetParInfo(current_par_ini);
		keys = std::get<0>(par_info);
		pars = std::get<1>(par_info);
		secc = std::get<2>(par_info);
		param_numb = std::get<3>(par_info);
		
		if (keys.size()>0) {
			show_pars = true;
		}
		bTryToLoad = false;
	}
	
	//Load ini if button is pressed
	if (ImGui::Button("Load par")) {
		//Loads the par into a INIStructure
		std::tuple<mINI::INIStructure, std::map<std::string, bool>> LoadedPar;
		LoadedPar = f_LoadPar(par_name);
		current_par_ini = std::get<0>(LoadedPar);
		aBool = std::get<1>(LoadedPar);
		//Sets keys, pars and secc
		std::tuple< std::vector<std::string>,std::vector<std::string>,std::vector<std::string>,std::vector<int> > par_info;
		par_info = GetParInfo(current_par_ini);
		keys = std::get<0>(par_info);
		pars = std::get<1>(par_info);
		secc = std::get<2>(par_info);
		param_numb = std::get<3>(par_info);
		
		show_pars = true;
	}
	
	//Show the parameters after clicking the "Load par" button
	if (show_pars) {
		/*
		i -> Moves between sections (secc)
		j -> Moves between keys (keys)
		k -> Moves between values (pars)
		-- keys.size() = pars.size()
		*/
		ImGui::InputText("##Output", par_name_to_export, 32);
		if (ImGui::Button("Save par")) {
			//Creates a new structure, saving the changes
			mINI::INIStructure new_par = CreatePar(secc, keys, pars, param_numb);
			//Generate the file
			f_SavePar(par_name_to_export, new_par);
		}
		int current_par = 0;
		//std::map<std::string, std::string> ParamAlias = std::get<0>(Aliases);
		//std::map<std::string, std::string> ToolTips = std::get<1>(Aliases);
		for (int i=0; i<secc.size(); i++) {
			ImGui::SeparatorText(secc[i].c_str());
			for (int j=0; j<param_numb[i]; j++) {
				if (ParamAlias.count(keys[current_par])) {
					ImGui::Text(ParamAlias[keys[current_par]].c_str()); //Displays the alias
				} else {
					ImGui::Text(keys[current_par].c_str()); //Displays the key
				}
				//ImGui::Text(keys[current_par].c_str()); //Displays the key
				ImGui::SameLine();
				std::string id_input = "##" + std::to_string(current_par);
				//USES DisplayParam from HERE CONSIDERING THE TYPES
				DisplayParam(keys[current_par], pars[current_par], id_input, current_par);
				current_par++;
			}
		}
		ImGui::SeparatorText("Execution");
		ImGui::Text("Executable name"); ImGui::SameLine(); HelpMarker("The compiled .exe file in the directory"); ImGui::SameLine();
		static std::string exe_name = GetExeName();
		static int np = std::thread::hardware_concurrency() - 1;
		if(np<1){np = 1;}
		ImGui::InputText("##ExeName", &exe_name);
		ImGui::Text("Number of CPU cores"); ImGui::SameLine(); HelpMarker("Number of cores where the simulation will run. By default uses almost all the cores that can detect, leaving one free");
		ImGui::SameLine();
		ImGui::InputInt("##NP", &np);
		if (ImGui::Button("Run simulation")) {
			//Runs MAD
			mINI::INIStructure new_par = CreatePar(secc, keys, pars, param_numb);
			f_RunMad(np, exe_name, new_par);
		}
		if (ImGui::Button("See results in PARAVIEW")) {
			//NEED A FIX. "paraview" can always run with this command???
			f_RunParaview("paraview", current_par_ini["IO"]["dirResults"], current_par_ini["IO"]["patientName"]);
		}
	}
	ImGui::End();
}

void display_GraphTest() {
	ImGui::Begin("Plots");
	
	ImGui::Text("Simulation results in real time. Hover the plots to see the values.");
	
	//static float values[90] = {};
	//static int values_offset = 0;
	//static double refresh_time = 0.0;
	//if (refresh_time == 0.0)
	//	refresh_time = ImGui::GetTime();
	//while (refresh_time < ImGui::GetTime()) // Create data at fixed 60 Hz rate for the demo
	//{
	//	values[values_offset] = f_GetDoubleFromBin("feedback.bin", 0); //Here goes the function value
	//	values_offset = (values_offset + 1) % IM_ARRAYSIZE(values);
	//	refresh_time += 1.0f / 60.0f;
	//}
	//ImGui::PlotLines("Incompressibility failure", values, IM_ARRAYSIZE(values), values_offset);
	
	ImPlot::CreateContext();
	int const size = 1000;
	std::vector<double> o_if;
	static double xs1[size], ys1[size];
    for (int i = 0; i < size; ++i) {
        xs1[i] = i; //* 0.001f;
        //ys1[i] = 0.5f + 0.5f * sinf(50 * (xs1[i] + (float)ImGui::GetTime() / 10));
		ys1[i] = f_GetDoubleFromBin("feedback.bin", 1);
    }
    
    if (ImPlot::BeginPlot("Imcompressibility failure")) {
        ImPlot::SetupAxes("t","I.F.");
		ImPlot::SetupAxesLimits(0,10,-0.001,0.001);
        ImPlot::PlotLine("Imcomp. Failure", xs1, ys1, size);
        ImPlot::EndPlot();
		//std::cout << "\n MMM current ys1: MMM" << ys1[0] << "\n";
    }
	ImPlot::DestroyContext();
	ImGui::End();
}

void display_DemoPlot() {
	ImPlot::CreateContext();
	ImPlot::ShowDemoWindow();
	ImPlot::DestroyContext();
}

class SubWindow {
	//Pendiente
public:
	std::string sWindowTitle;
	std::string sWelcomeText;
	
	void display();
};