/*=============================================================================
  This file is part of the code MAD
  Multi-physics for mechanicAl engineering and Data assimilation
  Copyright (C) 2017-2021,
    
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

#ifndef MAD_parameters
#define MAD_parameters

#include <STLvectorUtils.hpp>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <petsc.h>
#include <tools.hpp>
#include <regex>

using namespace std;

class Parameters {

  public:
    Parameters(){} /* Default constructor is needed to use parameters object as a member of other clases */
    Parameters(string dataFile);

    void print();
    void finalize();
    void finalize_short();

    /* Acces functions */
    inline const string & dissimilarity_matrix() const {
      return m_dissimilarity_matrix;}

    inline const string & errorMatrix() const {
      return m_errorMatrix;}

    inline const string & initial_condition_pre() const {
      return m_initial_condition_pre;}

    inline const string & initial_condition_vel() const {
      return m_initial_condition_vel;}

    inline const string & cartesian_grid() const {
      return m_cartesian_grid;}

    inline const string & voxelization_type() const {
      return m_voxelization_type;}

    inline const vector<int>& bdLabels() const {
      return m_bdLabels;}

    inline const vector<int>& bcNeumann() const {
      return m_bcNeumann;}

    inline const int & outlet() const {
      return m_outlet;}

    inline const int & jump() const {
      return m_jump;}

    inline const int & snaps_per_sim() const {
      return m_snaps_per_sim;}

    inline const int & noiseIterations() const {
      return m_noiseIterations;}

    inline const vector<int> & outlets() const {
      return m_outlets;}

    inline const int & nbVariables() const {
      return m_nbVariables;}

    inline const vector<int> & nbDofsPerNode() const {
      return m_nbDofsPerNode;}

    inline const vector<int> & measureIt() const {
      return m_measureIt;}

    inline const vector<double> & x_range() const {
      return m_x_range;}

    inline const vector<double> & factorInnerProduct() const {
      return m_factorInnerProduct;}

    inline const vector<double> & y_range() const {
      return m_y_range;}

    inline const vector<double> & z_range() const {
      return m_z_range;}

    inline const vector<string> & variableName() const {
      return m_variableName;}

    inline const vector<string> & fe() const {
      return m_fe;}

    inline const vector<double> & beamDir() const {
      return m_beamDir;}

    inline const vector<double> & transDir() const {
      return m_transDir;}

    inline const vector<double> & voxelSize() const {
      return m_voxelSize;}

    inline const int & nbVoxelsX() const {
      return m_nbVoxelsX;}

    inline const int & nbVoxelsY() const {
      return m_nbVoxelsY;}

    inline const int & nbVoxelsZ() const {
      return m_nbVoxelsZ;}

    inline const vector<int> & inlets() const {
      return m_inlets;}

    inline const vector<int> & walls() const {
      return m_walls;}

    inline const int & inlet() const {
      return m_inlet;}

    inline const int & wall() const {
      return m_wall;}

    inline const int & bdLabels(int i) const {
      return m_bdLabels[i];}

    inline const vector<double>& wave_velocity() const {
      return m_wave_velocity;}

    inline const vector<vector<double>>& bdConditions() const {
      return m_bdConditions;}

    inline const double & UStimeEnd() const {
      return m_UStimeEnd;}

    inline const double & scaleMesh() const {
      return m_scaleMesh;}

    inline const double & scaleUnits() const {
      return m_scaleUnits;}

    inline const double & termalExpansionCoefficient() const {
      return m_termalExpansionCoefficient;}

    inline const double & resistance() const {
      return m_resistance;}

    inline const double & damping() const {
      return m_damping;}

    inline const double & bingham() const {
      return m_bingham;}

    inline const double & amplitude() const {
      return m_amplitude;}

    inline const double & period() const {
      return m_period;}

    inline const double & ksp_tolerance() const {
      return m_ksp_tolerance;}

    inline const double & ksp_tolerance_div() const {
      return m_ksp_tolerance_div;}

    inline const int & ksp_max_iterations() const {
      return m_ksp_max_iterations;}

    inline const int & ksp_restartGMRESiterations() const {
      return m_ksp_restartGMRESiterations;}

    inline const double & ksp_tolerance_absolute() const {
      return m_ksp_tolerance_absolute;}

    inline const double & nonLinearTolerance() const {
      return m_nonLinearTolerance;}

    inline const double & modelCertainty() const {
      return m_modelCertainty;}

    inline const double & measuresCertainty() const {
      return m_measuresCertainty;}

    inline const vector<double> & resistances() const {
      return m_resistances;}

    inline const vector<double> & gravity() const {
      return m_gravity;}

    inline const double & frecuency() const {
      return m_frecuency;}

    // HAML added
    inline const double & backflowCoeff() const {
      return m_backflowCoeff;}

    inline const vector<double> & capacitances() const {
      return m_capacitances;}

    inline const vector<double> & distalResistances() const {
      return m_distalResistances;}

    inline const vector<double> & distalPressures0() const {
      return m_distalPressures0;}

    inline const double & inlet_radius() const {
      return m_inlet_radius;}

    inline const double & capacitance() const {
      return m_capacitance;}

    inline const double & distalResistance() const {
      return m_distalResistance;}

    inline const double & distalPress0() const {
      return m_distalPress0;}

    inline const double & inlet_u0() const {
      return m_inlet_u0;}

    inline const double & c() const {
      return m_c;}

    inline const double & diffusivity() const {
      return m_diffusivity;}

    inline const double & UStimeStart() const {
      return m_UStimeStart;}

    inline const double & UStimeStep() const {
      return m_UStimeStep;}

    inline const double & CFDtimeStep() const {
      return m_CFDtimeStep;}

    inline const int & nbModes() const {
      return m_nbModes;}

    inline const int & nbSamples() const {
      return m_nbSamples;}

    inline const string & dirMeasures() const {
      return m_dirMeasures;}

    inline const int & verbose() const {
      return m_verbose;}

    inline const string & inflow_data() const {
      return m_inflow_data;}

    inline const string & viscousTerm() const {
      return m_viscousTerm;}

    inline const string & solver() const {
      return m_solver;}

    inline const string & cfd_stab() const {
      return m_cfd_stab;}

    inline const bool & useModifiedGramSchmidt() const {
      return m_useModifiedGramSchmidt;}

    inline const string & preconditioner() const {
      return m_preconditioner;}

    inline const bool & convective() const {
      return m_convective;}

    inline const string & solutionFolder1() const {
      return m_solutionFolder1;}

    inline const string & solutionFolder2() const {
      return m_solutionFolder2;}

    inline const string & dirMatsStokes() const {
      return m_dirMatsStokes;}

    inline const string & surfaceMapping() const {
      return m_surfaceMapping;}

    inline const string & surfaceMappingInverse() const {
      return m_surfaceMappingInverse;}

    inline const string & dirMeasures2() const {
      return m_dirMeasures2;}

    inline const string & rieszRepresentersFolder() const {
      return m_rieszRepresentersFolder;}

    inline const string & geometryData() const {
      return m_geometryData;}

    inline const string & space_1() const {
      return m_space1;}

    inline const string & space_2() const {
      return m_space2;}

    inline const string & geometry() const {
      return m_geometry;}

    inline const string & dirModel() const {
      return m_dirModel;}

    inline const string & dirResults() const {
      return m_dirResults;}

    inline const string & measures() const {
      return m_measures;}

    inline const string & dirSyntheticField() const {
      return m_dirSyntheticField;}

    inline const double & gaussianNoiseLevel() const {
      return m_gaussianNoiseLevel;}

    inline const double & heartRate() const {
      return m_heartRate;}

    inline const string & modeMeasures() const {
      return m_modeMeasures;}

    inline const int & nbHeartRateWindows() const {
      return m_nbHeartRateWindows;}

    inline const int & nbTimeWindows() const {
      return m_nbTimeWindows;}

    inline const int & end() const {
      return m_end;}

    inline const int & start() const {
      return m_start;}

    inline const int & top() const {
      return m_top;}

    inline const int & bottom() const {
      return m_bottom;}

    inline const int & left() const {
      return m_left;}

    inline const int & right() const {
      return m_right;}

    inline const int & nonLinearMaxIterations() const {
      return m_nonLinearMaxIterations;}

    inline const int & nbVertices() const {
      return m_nbVertices;}

    inline const int & fixed_boundary() const {
      return m_fixed_boundary;}

    inline const size_t & nbMeasures() const {
      return m_nbMeasures;}

    inline const int & nbMeshes() const {
      return m_nbMeshes;}

    inline const string & outputFormat() const {
      return m_outputFormat;}

    inline const string & inputFormat() const {
      return m_inputFormat;}

    inline const string & measuresName() const {
      return m_measuresName;}

    inline const string & modelFormat() const {
      return m_modelFormat;}

    inline const bool & useSyntheticMeasures() const {
      return m_useSyntheticMeasures;}

    inline const bool & saveMatrices() const {
      return m_saveMatrices;}

    inline const bool & writeNonLinearIterations() const {
      return m_writeNonLinearIterations;}

    inline const bool & hemodynamics() const {
      return m_hemodynamics;}

    inline const bool & minimalOutput() const {
      return m_minimalOutput;}

    inline const bool & loadMeasuresAndG() const {
      return m_loadMeasuresAndG;}

    inline const bool & restrainedLeastSquares() const {
      return m_restrainedLeastSquares;}

    inline const bool & saveRieszRepresenters() const {
      return m_saveRieszRepresenters;}

    inline const bool & plotRieszRepresenters() const {
      return m_plotRieszRepresenters;}

    inline const bool & withMeasuresSimulator() const {
      return m_withMeasuresSimulator;}

    inline const double & USsampleVolumeSizeL() const {
      return m_USsampleVolumeSizeL;}

    inline const double & USsampleVolumeSizeT() const {
      return m_USsampleVolumeSizeT;}

    inline const double & USsampleVolumeSizeZ() const {
      return m_USsampleVolumeSizeZ;}

    inline const double & beamDirX() const {
      return m_beamDirX;}

    inline const double & beamDirY() const {
      return m_beamDirY;}

    inline const double & beamDirZ() const {
      return m_beamDirZ;}

    inline const double & transDirX() const {
      return m_transDirX;}

    inline const double & transDirY() const {
      return m_transDirY;}

    inline const double & transDirZ() const {
      return m_transDirZ;}

    inline const double & HRmax() const {
      return m_HRmax;}

    inline const double & HRmin() const {
      return m_HRmin;}

    inline const int & cardiacCycle() const {
      return m_cardiacCycle;}

    inline const string & patientName() const {
      return m_patientName;}

    inline const string & organ() const {
      return m_organ;}

    inline const string & method() const {
      return m_method;}

    inline const string & type() const {
      return m_type;}

    inline const string & innerProduct(int i = 0) const {
      return m_innerProduct[i];}

    inline const bool & plotTarget() const {
      return m_plotTarget;}

    inline const bool & convective_hammer() const {
      return m_convective_hammer;}

    inline const bool & monitorKSP() const {
      return m_monitorKSP;}

    inline const bool & mixed_problem() const {
      return m_mixed_problem;}

    inline const bool & backflowStab() const {
      return m_backflowStab;}

    inline const size_t & nbSimulations() const {
      return m_nbSimulations;}

    inline const size_t & nbQuadraturePoints() const {
      return m_nbQuadraturePoints;}

    inline const size_t & nbQuadraturePointsBD() const {
      return m_nbQuadraturePointsBD;}

    inline const size_t & nbQuadraturePointsTriangles() const {
      return m_nbQuadraturePointsTriangles;}

    inline const bool & plotMeasures() const {
      return m_plotMeasures;}

    inline const bool & reuse_preconditioner() const {
      return m_reuse_preconditioner;} 

    inline const bool & use_solution_as_guess_KSP() const {
      return m_use_solution_as_guess_KSP;}

    inline const bool & debug() const {
      return m_debug;}

    inline const bool & cleanLabels() const {
      return m_cleanLabels;}

    inline const bool & interruptor() const {
      return m_interruptor;}

    inline const bool & exportSurfaceOnly() const {
      return m_exportSurfaceOnly;}

    inline const bool & usePiola() const {
      return m_usePiola;}

    inline const bool & adaptive_nbModes() const {
      return m_adaptive_nbModes;}

    inline const bool & computeRieszRepresenters() const {
      return m_computeRieszRepresenters;}

    inline const string & templateGeometry() const {
      return m_templateGeometry;}

    inline const string & templateModel() const {
      return m_templateModel;}

    inline const string & model() const {
      return m_model;}

    inline const string & modeVectors() const {
      return m_modeVectors;}

    inline const string & interpolation_method() const {
      return m_interpolation_method;}

    inline const string & maniFolder() const {
      return m_maniFolder;}

    inline const string & stiffness_matrix_path() const {
      return m_stiffness_matrix_path;}

    inline const string & stiffness_matrix_path_template() const {
      return m_stiffness_matrix_path_template;}

    inline const string & mass_matrix_path() const {
      return m_mass_matrix_path;}

    inline const string & mass_scalar_matrix_path() const {
      return m_mass_scalar_matrix_path;}

    inline const string & mass_boundary_matrix_path() const {
      return m_mass_boundary_matrix_path;}

    inline const string & mass_matrix_path_pressure() const {
      return m_mass_matrix_path_pressure;}

    inline const string & problemType() const {
      return m_problemType;}

    inline const string & timeIntegration() const {
      return m_timeIntegration;}

    inline const double & timeStep() const {
      return m_timeStep;}

    inline const double & closureTime() const {
      return m_closureTime;}

    inline const size_t & nbIterations() const {
      return m_nbIterations;}

    inline const double & weight_stiffness() const {
      return m_weight_stiffness;}

    inline const double & weight_mass() const {
      return m_weight_mass;}

    inline const double & weight_pressure() const {
      return m_weight_pressure;}

    inline const double & viscosity() const {
      return m_viscosity;}

    inline const double & power_law_n() const {
      return m_power_law_n;}

    inline const double & power_law_m() const {
      return m_power_law_m;}

    inline const double & viscosity_solid() const {
      return m_viscosity_solid;}

    inline const double & viscosity_fluid() const {
      return m_viscosity_fluid;}

    inline const double & porosity() const {
      return m_porosity;}

    inline const double & fluid_density() const {
      return m_fluid_density;}

    inline const double & porosity2() const {
      return m_porosity2;}

    inline const double & density() const {
      return m_density;}

    inline const double & reynolds() const {
      return m_reynolds;}

    inline const double & inlet_symmetry() const {
      return m_inlet_symmetry;}

    inline const double & csf_pressure() const {
      return m_csf_pressure;}

    inline const double & csf_factor() const {
      return m_csf_factor;}

    inline const double & referenceTemperature() const {
      return m_referenceTemperature;}

    inline const double & alpha() const {
      return m_alpha;}

    inline const double & beta() const {
      return m_beta;}

    inline const double & gamma() const {
      return m_gamma;}

    inline const double & delta() const {
      return m_delta;}

    inline const double & epsilon() const {
      return m_epsilon;}

    inline const double & zeta() const {
      return m_zeta;}

    inline const double & eta() const {
      return m_eta;}

    inline const double & theta() const {
      return m_theta;}

    inline const double & iota() const {
      return m_iota;}

    inline const double & kappa() const {
      return m_kappa;}

    inline const double & lambda() const {
      return m_lambda;}

    inline const double & mu() const {
      return m_mu;}

    inline const double & nu() const {
      return m_nu;}

    inline const double & omicron() const {
      return m_omicron;}

    inline const double & pi() const {
      return m_pi;}

    inline const double & rho() const {
      return m_rho;}

    inline const double & sigma() const {
      return m_sigma;}

    inline const double & tau() const {
      return m_tau;}

    inline const double & upsilon() const {
      return m_upsilon;}

    inline const double & phi() const {
      return m_phi;}

    inline const double & chi() const {
      return m_chi;}

    inline const double & tsi() const {
      return m_tsi;}

    inline const double & omega() const {
      return m_omega;}

    inline const vector<double> & poissonRatio() const {
      return m_poissonRatio;}

    inline const vector<double> & youngModulus() const {
      return m_youngModulus;}

    inline const vector<double> & permeability() const {
      return m_permeability;}

    inline const double & heigth() const {
      return m_heigth;}

    inline const double & length() const {
      return m_length;}

    inline const double & diameter() const {
      return m_diameter;}

  private:

    void assertParameters();

    bool m_debug = false;

    /* Time */
    size_t m_nbIterations;
    int m_start;
    int m_end;

    /* Finite elements */
    size_t m_nbQuadraturePoints;
    size_t m_nbQuadraturePointsBD;
    size_t m_nbQuadraturePointsTriangles;
    bool m_mixed_problem = false;
    vector<int> m_bdLabels;
    vector<int> m_bcNeumann;
    vector<vector<double>> m_bdConditions;
    int m_nbVariables = 1;
    vector<int> m_nbDofsPerNode;
    vector<string> m_variableName;
    vector<string> m_fe;

    /* Hilbert spaces */
    string m_problemType = "vector";
    vector<string> m_innerProduct;
    string m_mass_matrix_path = "";
    string m_mass_scalar_matrix_path = "";
    string m_mass_boundary_matrix_path = "";
    string m_mass_matrix_path_pressure = "";
    string m_stiffness_matrix_path = "";
    double m_weight_stiffness = 1.0;
    double m_weight_mass = 1.0;
    double m_weight_pressure = 1.0;
    bool m_saveMatrices = false;
    bool m_hemodynamics = false;
    string m_inflow_data = "";

    /* Model */
    int m_nbTimeWindows = 1; // linear by default
    int m_nbHeartRateWindows = 1; // linear by default
    string m_dirModel;
    string m_modelFormat = ".vct";
    string m_model;
    string m_templateModel;
    int m_nbModes;
    double m_HRmax, m_HRmin;
    string m_method;
    string m_type = "linear";
    int m_snaps_per_sim = -1;
    int m_jump = 1;

    /* Model (offline) */
    string m_maniFolder;
    double m_timeStep;
    size_t m_nbSimulations;
    double m_heartRate;
    int m_cardiacCycle;
    string m_patientName;
    string m_organ;
    string m_measures;

    /* Geometry */
    int m_nbVertices;
    string m_geometryData;
    string m_surfaceMapping;
    string m_surfaceMappingInverse;
    string m_geometry = "";
    bool m_usePiola;
    bool m_exportSurfaceOnly = false;
    bool m_cleanLabels = false;
    string m_templateGeometry;
    string m_stiffness_matrix_path_template;
    string m_interpolation_method;
    int m_nbMeshes;
    int m_fixed_boundary;
    string m_space1;
    string m_space2;
    double m_scaleMesh = 1.0;
    double m_scaleUnits = 1.0;
    int m_top, m_bottom, m_left, m_right;
    
    /* Ultrasound */
    string m_dirMeasures;
    string m_dirMeasures2;
    string m_measuresName;
    string m_modeMeasures;
    string m_cfd_stab;
    double m_UStimeStart;
    double m_UStimeStep;
    double m_CFDtimeStep;
    double m_UStimeEnd;
    double m_USsampleVolumeSizeL;
    double m_USsampleVolumeSizeT;
    double m_USsampleVolumeSizeZ;
    vector<double> m_beamDir;
    vector<double> m_factorInnerProduct;
    vector<double> m_transDir;
    vector<double> m_voxelSize;
    int m_nbVoxelsX;
    int m_nbVoxelsY;
    int m_nbVoxelsZ;
    double m_beamDirX;
    double m_beamDirY;
    double m_beamDirZ;
    double m_transDirX;
    double m_transDirY;
    double m_transDirZ;
    int m_verbose = 5;
    size_t m_nbMeasures;

    /* Learning */
    string m_voxelization_type; 
    string m_cartesian_grid; 
    string m_dissimilarity_matrix;
    string m_errorMatrix;
    int m_noiseIterations;

    /* Ultrasound synthetic */
    bool m_useSyntheticMeasures = true;
    string m_dirSyntheticField;
    double m_gaussianNoiseLevel;
    bool m_withMeasuresSimulator;
    vector<int> m_measureIt;

    /* Input/Output */
    string m_dirResults;
    string m_outputFormat = ".vct"; 
    string m_inputFormat = ".bin";
    bool m_plotRieszRepresenters = false; // this saves .vct for visualization
    bool m_saveRieszRepresenters = false; // this saves .bin for future computations
    bool m_plotMeasures = false;
    bool m_plotTarget = false;
    bool m_writeNonLinearIterations = false;
    bool m_loadMeasuresAndG = false;
    bool m_computeRieszRepresenters = true;
    string m_rieszRepresentersFolder = " ";

    /* CFD */
    string m_viscousTerm = "symmetric";
    vector<double> m_gravity;
    bool m_backflowStab;
    double m_backflowCoeff = 0.2; // HAML added
    double m_viscosity;
    double m_viscosity_solid;
    double m_viscosity_fluid;
    double m_density;
    double m_reynolds;
    double m_inlet_symmetry;
    double m_inlet_u0 = 1.0;
    string m_dirMatsStokes;
    int m_inlet, m_outlet, m_wall;
    vector<int> m_inlets, m_walls, m_outlets;
    double m_inlet_radius;
    double m_power_law_m = 1.0;
    double m_power_law_n = 1.0;
    double m_termalExpansionCoefficient;
    double m_referenceTemperature;
    bool m_convective = true;

    /* greek */
    double m_alpha;
    double m_beta;
    double m_gamma;
    double m_delta;
    double m_epsilon;
    double m_zeta;
    double m_eta;
    double m_theta;
    double m_iota;
    double m_kappa;
    double m_lambda;
    double m_mu;
    double m_nu;
    double m_xi;
    double m_omicron;
    double m_pi;
    double m_rho;
    double m_sigma;
    double m_tau;
    double m_upsilon;
    double m_phi;
    double m_chi;
    double m_tsi;
    double m_omega;

    /* poro-elasticity */
    double m_porosity;
    vector<double> m_permeability;
    double m_fluid_density;
    double m_porosity2;
    vector<double> m_youngModulus;
    vector<double> m_poissonRatio;
    double m_period;
    double m_amplitude;
    double m_frecuency;
    double m_csf_pressure;
    double m_csf_factor;
    double m_diffusivity;
    double m_closureTime;
    double m_c;

    /* Windkessel */
    double m_resistance;
    double m_capacitance;
    double m_distalResistance;
    double m_distalPress0;
    vector<double> m_resistances;
    vector<double> m_capacitances;
    vector<double> m_distalResistances;
    vector<double> m_distalPressures0;

    bool m_restrainedLeastSquares;
    bool m_interruptor;
    bool m_adaptive_nbModes;
 
    int m_nbSamples;
    double m_modelCertainty;
    double m_measuresCertainty;

    /* Solver */
    string m_initial_condition_pre = "none";
    string m_initial_condition_vel = "none";
    string m_solver;
    bool  m_useModifiedGramSchmidt = true;
    string m_preconditioner;
    string m_timeIntegration = "BDF1";
    bool m_use_solution_as_guess_KSP = false;
    bool m_reuse_preconditioner = false;
    int m_ksp_max_iterations = 10000;
    int m_ksp_restartGMRESiterations = 10000;
    double m_ksp_tolerance = 1e-2;
    double m_ksp_tolerance_absolute = 1e-2;
    double m_ksp_tolerance_div = 1e-2;
    double m_nonLinearTolerance = 999; /* Linear by default */
    int m_nonLinearMaxIterations = 1;  /* Linear by default */
    bool m_monitorKSP;

    /* Hammer */
    double m_bingham;
    double m_diameter;
    double m_length;
    double m_heigth;
    bool m_convective_hammer;
   
    /* folding */ 
    string m_solutionFolder1;
    string m_solutionFolder2;
    string m_name;
    bool m_minimalOutput = false;
    string m_modeVectors = "contiguous";

    void deploy_vector(vector<int> & vector_to_deploy, string parameter_value);
    void deploy_vector(vector<double> & vector_to_deploy, string parameter_value);
    void deploy_vector(vector<string> & vector_to_deploy, string parameter_value);
    void deploy_matrix(vector<vector<double>> & matrix_to_deploy, string parameter_value);

    double m_damping;
    vector<double> m_wave_velocity;
    vector<double> m_x_range;
    vector<double> m_y_range;
    vector<double> m_z_range;

    int m_world_rank;  

    Timer reloj;

};  

#endif
