/*=============================================================================
  This file is part of the code MAD
  Multi-physics for mechanicAl engineering and Data assimilation
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

#include <parameters.hpp>

using namespace std;

Parameters::Parameters(string dataFile){

  m_name = dataFile;

  ifstream parametersFile(dataFile);

  /* Set some default variables so the code don't crash. It is better to assert and exit? */
  m_variableName.push_back("no_named_variable");

  if (parametersFile.is_open()){

    string line;
    
    while (getline(parametersFile,line)){ 

      string parameter_name = line.substr(0, line.find("="));
      string parameter_value = line.substr(line.find("=")+1, line.length());

      /* Replacing os enviroment variables */
      parameter_value = std::regex_replace(parameter_value, std::regex("MAD_DATA"), get_env_var("MAD_DATA"));
      parameter_value = std::regex_replace(parameter_value, std::regex("MAD_RESULTS"), get_env_var("MAD_RESULTS"));

      if (parameter_name == "nbTimeWindows")
        m_nbTimeWindows = stoi(parameter_value);

      if (parameter_name == "nbVariables")
        m_nbVariables = stoi(parameter_value);
      
      if (parameter_name == "nbHeartRateWindows")
        m_nbHeartRateWindows = stoi(parameter_value);

      if (parameter_name == "templateGeometry")
        m_templateGeometry = parameter_value;

      if (parameter_name == "dissimilarity_matrix")
        m_dissimilarity_matrix = parameter_value;

      if (parameter_name == "errorMatrix")
        m_errorMatrix = parameter_value;

      if (parameter_name == "interpolation_method")
        m_interpolation_method = parameter_value;

      if (parameter_name == "model")
        m_model = parameter_value;

      if (parameter_name == "cfd_stab")
        m_cfd_stab = parameter_value;

      if (parameter_name == "templateModel")
        m_templateModel = parameter_value;

      if (parameter_name == "inflow_data")
        m_inflow_data = parameter_value;

      if (parameter_name == "dirModel")
        m_dirModel = parameter_value;

      if (parameter_name == "surfaceMapping")
        m_surfaceMapping = parameter_value;

      if (parameter_name == "surfaceMappingInverse")
        m_surfaceMappingInverse = parameter_value;

      if (parameter_name == "cardiacCycle")
        m_cardiacCycle = stoi(parameter_value);

      if (parameter_name == "noiseIterations")
        m_noiseIterations = stoi(parameter_value);

      if (parameter_name == "start")
        m_start = stoi(parameter_value);

      if (parameter_name == "nonLinearMaxIterations")
        m_nonLinearMaxIterations = stoi(parameter_value);

      if (parameter_name == "ksp_restartGMRESiterations")
        m_ksp_restartGMRESiterations= stoi(parameter_value);

      if (parameter_name == "nbVoxelsX")
        m_nbVoxelsX = stoi(parameter_value);

      if (parameter_name == "nbVoxelsY")
        m_nbVoxelsY = stoi(parameter_value);

      if (parameter_name == "nbVoxelsZ")
        m_nbVoxelsZ = stoi(parameter_value);

      if (parameter_name == "end")
        m_end = stoi(parameter_value);

      if (parameter_name == "nbSamples")
        m_nbSamples = stoi(parameter_value);

      if (parameter_name == "nbModes")
        m_nbModes = stoi(parameter_value);

      if (parameter_name == "density")
        m_density = stod(parameter_value);

      if (parameter_name == "damping")
        m_damping = stod(parameter_value);

      if (parameter_name == "bingham")
        m_bingham = stod(parameter_value);

      if (parameter_name == "reynolds")
        m_reynolds = stod(parameter_value);

      if (parameter_name == "csf_factor")
        m_csf_factor = stod(parameter_value);

      if (parameter_name == "inlet_radius")
        m_inlet_radius = stod(parameter_value);

      if (parameter_name == "inlet_u0")
        m_inlet_u0 = stod(parameter_value);

      if (parameter_name == "amplitude")
        m_amplitude = stod(parameter_value);

      // HAML added: beta for Neumann backflow stabilization
      if (parameter_name == "backflowCoeff")
        m_backflowCoeff = stod(parameter_value);

      if (parameter_name == "period")
        m_period = stod(parameter_value);

      if (parameter_name == "modelCertainty")
        m_modelCertainty = stod(parameter_value);

      if (parameter_name == "measuresCertainty")
        m_measuresCertainty = stod(parameter_value);

      if (parameter_name == "inlet_symmetry")
        m_inlet_symmetry = stod(parameter_value);

      if (parameter_name == "viscosity")
        m_viscosity = stod(parameter_value);

      if (parameter_name == "power_law_n")
        m_power_law_n = stod(parameter_value);

      if (parameter_name == "power_law_m")
        m_power_law_m = stod(parameter_value);

      if (parameter_name == "viscosity_fluid")
        m_viscosity_fluid = stod(parameter_value);

      if (parameter_name == "viscosity_solid")
        m_viscosity_solid = stod(parameter_value);

      if (parameter_name == "porosity")
        m_porosity = stod(parameter_value);

      if (parameter_name == "length")
        m_length = stod(parameter_value);

      if (parameter_name == "diameter")
        m_diameter = stod(parameter_value);

      if (parameter_name == "heigth")
        m_heigth = stod(parameter_value);

      if (parameter_name == "permeability")
        deploy_vector(m_permeability, parameter_value);

      if (parameter_name == "fluid_density")
        m_fluid_density = stod(parameter_value);

      if (parameter_name == "porosity2")
        m_porosity2 = stod(parameter_value);

      if (parameter_name == "c")
        m_c = stod(parameter_value);

      if (parameter_name == "diffusivity")
        m_diffusivity = stod(parameter_value);

      if (parameter_name == "heartRate")
        m_heartRate = stod(parameter_value);

      if (parameter_name == "nbVertices")
        m_nbVertices = stoi(parameter_value);

      if (parameter_name == "nbIterations")
        m_nbIterations = stoi(parameter_value);

      if (parameter_name == "nbMeasures")
        m_nbMeasures = stoi(parameter_value);

      if (parameter_name == "nbMeshes")
        m_nbMeshes = stoi(parameter_value);

      if (parameter_name == "nbQuadraturePoints")
        m_nbQuadraturePoints = stoi(parameter_value);

      if (parameter_name == "nbQuadraturePointsBD")
        m_nbQuadraturePointsBD = stoi(parameter_value);

      if (parameter_name == "nbQuadraturePointsTriangles")
        m_nbQuadraturePointsTriangles = stoi(parameter_value);

      if (parameter_name == "geometryData"){
        m_geometryData = parameter_value;
      }

      if (parameter_name == "modeVectors"){
        m_modeVectors = parameter_value;
      }

      if (parameter_name == "space_1")
        m_space1 = parameter_value;

      if (parameter_name == "space_2")
        m_space2 = parameter_value;

      if (parameter_name == "geometry")
        m_geometry = parameter_value;

      if (parameter_name == "measuresName")
        m_measuresName = parameter_value;

      if (parameter_name == "dirMeasures")
        m_dirMeasures = parameter_value;

      if (parameter_name == "top")
        m_top = stoi(parameter_value);

      if (parameter_name == "bottom")
        m_bottom = stoi(parameter_value);

      if (parameter_name == "left")
        m_left = stoi(parameter_value);

      if (parameter_name == "right")
        m_right = stoi(parameter_value);

      if (parameter_name == "verbose")
        m_verbose = stoi(parameter_value);

      if (parameter_name == "fixed_boundary")
        m_fixed_boundary = stoi(parameter_value);

      if (parameter_name == "preconditioner")
        m_preconditioner = parameter_value;

      if (parameter_name == "solutionFolder1")
        m_solutionFolder1 = parameter_value;

      if (parameter_name == "solutionFolder2")
        m_solutionFolder2 = parameter_value;

      if (parameter_name == "solver")
        m_solver = parameter_value;

      if (parameter_name == "timeIntegration")
        m_timeIntegration = parameter_value;

      if (parameter_name == "viscousTerm")
        m_viscousTerm = parameter_value;

      if (parameter_name == "closureTime")
        m_closureTime = stod(parameter_value);

      if (parameter_name == "dirMatsStokes")
        m_dirMatsStokes = parameter_value;

      if (parameter_name == "dirMeasures2")
        m_dirMeasures2 = parameter_value;

      if (parameter_name == "rieszRepresentersFolder")
        m_rieszRepresentersFolder = parameter_value;

      if (parameter_name == "beamDirX")
        m_beamDirX = stod(parameter_value);

      if (parameter_name == "beamDirY")
        m_beamDirY = stod(parameter_value);

      if (parameter_name == "beamDirZ")
        m_beamDirZ = stod(parameter_value);

      if (parameter_name == "transDirX")
        m_transDirX = stod(parameter_value);

      if (parameter_name == "transDirY")
        m_transDirY = stod(parameter_value);

      if (parameter_name == "transDirZ")
        m_transDirZ = stod(parameter_value);

      if (parameter_name == "modeMeasures")
        m_modeMeasures = parameter_value;

      if (parameter_name == "UStimeStart")
        m_UStimeStart = stod(parameter_value);

      if (parameter_name == "USsampleVolumeSizeL")
        m_USsampleVolumeSizeL = stod(parameter_value);

      if (parameter_name == "USsampleVolumeSizeT")
        m_USsampleVolumeSizeT = stod(parameter_value);

      if (parameter_name == "USsampleVolumeSizeZ")
        m_USsampleVolumeSizeZ = stod(parameter_value);

      if (parameter_name == "UStimeEnd")
        m_UStimeEnd = stod(parameter_value);

      if (parameter_name == "UStimeStep")
        m_UStimeStep = stod(parameter_value);

      if (parameter_name == "CFDtimeStep")
        m_CFDtimeStep = stod(parameter_value);

      if (parameter_name == "plotRieszRepresenters"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_plotRieszRepresenters = true;
        } else {
          m_plotRieszRepresenters = false;
        }
      }

      if (parameter_name == "minimalOutput"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_minimalOutput = true;
        } else {
          m_minimalOutput = false;
        }
      }

      if (parameter_name == "restrainedLeastSquares"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_restrainedLeastSquares = true;
        } else {
          m_restrainedLeastSquares = false;
        }
      }

      if (parameter_name == "adaptive_nbModes"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_adaptive_nbModes = true;
        } else {
          m_adaptive_nbModes = false;
        }
      }

      if (parameter_name == "loadMeasuresAndG"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_loadMeasuresAndG = true;
        } else {
          m_loadMeasuresAndG = false;
        }
      }

      if (parameter_name == "saveMatrices"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_saveMatrices = true;
        } else {
          m_saveMatrices = false;
        }
      }

      if (parameter_name == "hemodynamics"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_hemodynamics = true;
        } else {
          m_hemodynamics = false;
        }
      }

      if (parameter_name == "usePiola"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_usePiola = true;
        } else {
          m_usePiola = false;
        }
      }

      if (parameter_name == "convective"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_convective = true;
        } else {
          m_convective = false;
        }
      }

      if (parameter_name == "useModifiedGramSchmidt"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_useModifiedGramSchmidt = true;
        } else {
          m_useModifiedGramSchmidt = false;
        }
      }


      if (parameter_name == "use_solution_as_guess_KSP"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_use_solution_as_guess_KSP = true;
        } else {
          m_use_solution_as_guess_KSP = false;
        }
      }

      if (parameter_name == "reuse_preconditioner"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_reuse_preconditioner = true;
        } else {
          m_reuse_preconditioner = false;
        }
      }

      if (parameter_name == "writeNonLinearIterations"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_writeNonLinearIterations = true;
        } else {
          m_writeNonLinearIterations = false;
        }
      }

      if (parameter_name == "interruptor"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_interruptor = true;
        } else {
          m_interruptor = false;
        }
      }

      if (parameter_name == "exportSurfaceOnly"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_exportSurfaceOnly = true;
        } else {
          m_exportSurfaceOnly = false;
        }
      }

      if (parameter_name == "plotTarget"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_plotTarget = true;
        } else {
          m_plotTarget = false;
        }
      }

      if (parameter_name == "convective_hammer"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_convective_hammer = true;
        } else {
          m_convective_hammer = false;
        }
      }

      if (parameter_name == "debug"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_debug = true;
        } else {
          m_debug = false;
        }
      }

      if (parameter_name == "cleanLabels"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_cleanLabels = true;
        } else {
          m_cleanLabels = false;
        }
      }

      if (parameter_name == "monitorKSP"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_monitorKSP = true;
        } else {
          m_monitorKSP = false;
        }
      }

      if (parameter_name == "backflowStab"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_backflowStab = true;
        } else {
          m_backflowStab = false;
        }
      }

      if (parameter_name == "mixed_problem"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_mixed_problem = true;
        } else {
          m_mixed_problem = false;
        }
      }

      if (parameter_name == "computeRieszRepresenters"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_computeRieszRepresenters = true;
        } else {
          m_computeRieszRepresenters = false;
        }
      }

      if (parameter_name == "plotMeasures"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_plotMeasures = true;
        } else {
          m_plotMeasures = false;
        }
      }

      if (parameter_name == "useSyntheticMeasures"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_useSyntheticMeasures = true;
        } else {
          m_useSyntheticMeasures = false;
        }
      }

      if (parameter_name == "saveRieszRepresenters"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_saveRieszRepresenters = true;
        } else {
          m_saveRieszRepresenters = false;
        }
      }

      if (parameter_name == "withMeasuresSimulator"){
        if (parameter_value == "true" || parameter_value == "True" || parameter_value == "1" || parameter_value == "yes" || parameter_value == "Yes" || parameter_value == "YES"){
          m_withMeasuresSimulator = true;
        } else {
          m_withMeasuresSimulator = false;
        }
      }

      if (parameter_name == "dirSyntheticField")
        m_dirSyntheticField = parameter_value;

      if (parameter_name == "gaussianNoiseLevel")
        m_gaussianNoiseLevel = stod(parameter_value);

      if (parameter_name == "HRmax")
        m_HRmax = stod(parameter_value);

      if (parameter_name == "HRmin")
        m_HRmin = stod(parameter_value);

      if (parameter_name == "ksp_tolerance")
        m_ksp_tolerance = stod(parameter_value);

      if (parameter_name == "ksp_tolerance_div")
        m_ksp_tolerance_div = stod(parameter_value);

      if (parameter_name == "ksp_max_iterations")
        m_ksp_max_iterations = stod(parameter_value);

      if (parameter_name == "ksp_tolerance_absolute")
        m_ksp_tolerance_absolute = stod(parameter_value);

      if (parameter_name == "nonLinearTolerance")
        m_nonLinearTolerance = stod(parameter_value);

      if (parameter_name == "dirResults"){
        m_dirResults = parameter_value;
      }

      if (parameter_name == "initial_condition_vel")
        m_initial_condition_vel = parameter_value;

      if (parameter_name == "initial_condition_pre")
        m_initial_condition_pre = parameter_value;

      if (parameter_name == "measures")
        m_measures = parameter_value;

      if (parameter_name == "patientName")
        m_patientName = parameter_value;

      if (parameter_name == "organ")
        m_organ = parameter_value;

      if (parameter_name == "outputFormat")
        m_outputFormat = parameter_value;

      if (parameter_name == "inputFormat")
        m_inputFormat = parameter_value;

      if (parameter_name == "modelFormat")
        m_modelFormat = parameter_value;

      if (parameter_name == "method")
        m_method = parameter_value;

      if (parameter_name == "type")
        m_type = parameter_value;

      if (parameter_name == "maniFolder")
        m_maniFolder = parameter_value;

      if (parameter_name == "snaps_per_sim")
        m_snaps_per_sim = stoi(parameter_value);

      if (parameter_name == "innerProduct")
        deploy_vector(m_innerProduct, parameter_value);

      if (parameter_name == "measureIt")
        deploy_vector(m_measureIt, parameter_value);

      if (parameter_name == "weight_stiffness")
        m_weight_stiffness = stod(parameter_value);

      if (parameter_name == "weight_mass")
        m_weight_mass = stod(parameter_value);

      if (parameter_name == "csf_pressure")
        m_csf_pressure = stod(parameter_value);

      if (parameter_name == "termalExpansionCoefficient")
        m_termalExpansionCoefficient = stod(parameter_value);

      if (parameter_name == "weight_pressure")
        m_weight_pressure = stod(parameter_value);

      if (parameter_name == "distalResistance")
        m_distalResistance = stod(parameter_value);

      if (parameter_name == "distalPress0")
        m_distalPress0 = stod(parameter_value);

      if (parameter_name == "capacitance")
        m_capacitance = stod(parameter_value);

      if (parameter_name == "alpha")
        m_alpha = stod(parameter_value);

      if (parameter_name == "beta")
        m_beta = stod(parameter_value);

      if (parameter_name == "gamma")
        m_gamma = stod(parameter_value);

      if (parameter_name == "delta")
        m_delta = stod(parameter_value);

      if (parameter_name == "epsilon")
        m_epsilon = stod(parameter_value);

      if (parameter_name == "zeta")
        m_zeta = stod(parameter_value);

      if (parameter_name == "eta")
        m_eta = stod(parameter_value);

      if (parameter_name == "theta")
        m_theta = stod(parameter_value);

      if (parameter_name == "iota")
        m_iota = stod(parameter_value);

      if (parameter_name == "kappa")
        m_kappa = stod(parameter_value);

      if (parameter_name == "lambda")
        m_lambda = stod(parameter_value);

      if (parameter_name == "mu")
        m_mu = stod(parameter_value);

      if (parameter_name == "nu")
        m_nu = stod(parameter_value);

      if (parameter_name == "xi")
        m_xi = stod(parameter_value);

      if (parameter_name == "omicron")
        m_omicron = stod(parameter_value);

      if (parameter_name == "pi")
        m_pi = stod(parameter_value);

      if (parameter_name == "phi")
        m_phi = stod(parameter_value);

      if (parameter_name == "rho")
        m_rho = stod(parameter_value);

      if (parameter_name == "sigma")
        m_sigma = stod(parameter_value);

      if (parameter_name == "tau")
        m_tau = stod(parameter_value);

      if (parameter_name == "upsilon")
        m_upsilon = stod(parameter_value);

      if (parameter_name == "phi")
        m_phi = stod(parameter_value);

      if (parameter_name == "chi")
        m_chi = stod(parameter_value);

      if (parameter_name == "tsi")
        m_tsi = stod(parameter_value);

      if (parameter_name == "omega")
        m_omega = stod(parameter_value);

      if (parameter_name == "poissonRatio")
        deploy_vector(m_poissonRatio, parameter_value);

      if (parameter_name == "gravity")
        deploy_vector(m_gravity, parameter_value);

      if (parameter_name == "frecuency")
        m_frecuency = stod(parameter_value);

      if (parameter_name == "scaleMesh")
        m_scaleMesh = stod(parameter_value);

      if (parameter_name == "scaleUnits")
        m_scaleUnits = stod(parameter_value);

      if (parameter_name == "referenceTemperature")
        m_referenceTemperature = stod(parameter_value);

      if (parameter_name == "resistance")
        m_resistance = stod(parameter_value);

      if (parameter_name == "mass_matrix_path")
        m_mass_matrix_path = parameter_value;

      if (parameter_name == "mass_scalar_matrix_path")
        m_mass_scalar_matrix_path = parameter_value;

      if (parameter_name == "mass_boundary_matrix_path")
        m_mass_boundary_matrix_path = parameter_value;

      if (parameter_name == "mass_matrix_path_pressure")
        m_mass_matrix_path_pressure = parameter_value;

      if (parameter_name == "problemType")
        m_problemType = parameter_value;

      if (parameter_name == "cartesian_grid")
        m_cartesian_grid = parameter_value;

      if (parameter_name == "voxelization_type")
        m_voxelization_type = parameter_value;

      if (parameter_name == "stiffness_matrix_path_template")
        m_stiffness_matrix_path_template = parameter_value;

      if (parameter_name == "stiffness_matrix_path")
        m_stiffness_matrix_path = parameter_value;

      if (parameter_name == "timeStep")
        m_timeStep = stod(parameter_value);

      if (parameter_name == "nbSimulations")
        m_nbSimulations = stoi(parameter_value);

      if (parameter_name == "jump")
        m_jump = stoi(parameter_value);

      if (parameter_name == "inlet")
        m_inlet = stoi(parameter_value);

      if (parameter_name == "outlet")
        m_outlet = stoi(parameter_value);

      if (parameter_name == "outlets")
        deploy_vector(m_outlets, parameter_value);

      if (parameter_name == "factorInnerProduct")
        deploy_vector(m_factorInnerProduct, parameter_value);

      if (parameter_name == "youngModulus")
        deploy_vector(m_youngModulus, parameter_value);

      if (parameter_name == "bcNeumann")
        deploy_vector(m_bcNeumann, parameter_value);

      if (parameter_name == "x_range")
        deploy_vector(m_x_range, parameter_value);

      if (parameter_name == "y_range")
        deploy_vector(m_y_range, parameter_value);

      if (parameter_name == "z_range")
        deploy_vector(m_z_range, parameter_value);

      if (parameter_name == "beamDir")
        deploy_vector(m_beamDir, parameter_value);

      if (parameter_name == "voxelSize")
        deploy_vector(m_voxelSize, parameter_value);

      if (parameter_name == "transDir")
        deploy_vector(m_transDir, parameter_value);

      if (parameter_name == "inlets")
        deploy_vector(m_inlets, parameter_value);

      if (parameter_name == "walls")
        deploy_vector(m_walls, parameter_value);

      if (parameter_name == "distalResistances")
        deploy_vector(m_distalResistances, parameter_value);

      if (parameter_name == "resistances")
        deploy_vector(m_resistances, parameter_value);

      if (parameter_name == "capacitances")
        deploy_vector(m_capacitances, parameter_value);

      if (parameter_name == "distalPressures0")
        deploy_vector(m_distalPressures0, parameter_value);

      if (parameter_name == "wall")
        m_wall = stoi(parameter_value);

      if (parameter_name == "bdLabels")
        deploy_vector(m_bdLabels, parameter_value);

      if (parameter_name == "nbDofsPerNode")
        deploy_vector(m_nbDofsPerNode, parameter_value);

      if (parameter_name == "wave_velocity")
        deploy_vector(m_wave_velocity, parameter_value);

      if (parameter_name == "variableName")
        deploy_vector(m_variableName, parameter_value);

      if (parameter_name == "fe")
        deploy_vector(m_fe, parameter_value);

      if (parameter_name == "bdConditions")
        deploy_matrix(m_bdConditions, parameter_value);

    }

  } else {
    cout << "ERROR: can't open parameters file " << dataFile << endl; 
    exit(1);
  }

  /* by default, inner product matrix is not weigthed */
  if (m_innerProduct.size() > 0 & m_factorInnerProduct.size() == 0){
    m_factorInnerProduct.resize(m_innerProduct.size(), 1.0);
  }

  parametersFile.close();

//  assertParameters();

}


void Parameters::assertParameters(){
  assert(m_bdLabels.size() == m_bdConditions.size() && "Boundary conditions and labels should have the same size.");
  assert(m_problemType == "vector" or m_problemType == "scalar");
}

void Parameters::print(){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank != 0){
    return;
  }
  cout << endl;
  cout << "\n - - - - - - - - - "<< "\033[1;31mM.A.D.\033[0m" << " - - - - - - - - - " << endl;
  cout << endl;
  cout << "\tFelipe Galarce. 2017/2024." << endl;
//  cout << "\t" << show_date_();
  cout << "\n - - - - - - - - - - - - - - - - - - - - - - " << endl;
  cat(m_name);
  cout << "\n - - - - - - - - - - - - - - - - - - - - - - " << endl << endl;

}

void Parameters::finalize(){
  if (m_world_rank != 0){
    return;
  }
  cout << "\n - - - - - - - - - - Normal end - - - - - - - - - -" << endl;
  time_t now = time(0);
  cout << "\t " << ctime(&now);
  cout << "\t Total execution time: " << reloj.toc() << " seconds " << endl;
  cout << "\t Check and post-process results using: " << endl;
  cout << "\t paraview " << m_dirResults << "/" << m_patientName << ".case &" << endl;
  cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
}

void Parameters::finalize_short(){
  if (m_world_rank != 0){
    return;
  }
  cout << "\n - - - - - - - - - - Normal end - - - - - - - - - -" << endl;
  time_t now = time(0);
  cout << "\t " << ctime(&now);
  cout << "\t Total execution time: " << reloj.toc() << " seconds " << endl;
  cout << " - - - - - - - - - - - - - - - - - - - - - - - - - - " << endl;
}

void Parameters::deploy_vector(vector<int> & vector_to_deploy, string parameter_value){
  vector_to_deploy.clear();
  while (parameter_value.find(" ") != string::npos){ 
    vector_to_deploy.push_back( stoi( parameter_value.substr(0, parameter_value.find(" ")) ) );
    parameter_value = parameter_value.substr(parameter_value.find(" ")+1, parameter_value.length());
  }
  vector_to_deploy.push_back(stoi( parameter_value ));
}

void Parameters::deploy_vector(vector<double> & vector_to_deploy, string parameter_value){
  vector_to_deploy.clear();
  while (parameter_value.find(" ") != string::npos){ 
    vector_to_deploy.push_back( stod( parameter_value.substr(0, parameter_value.find(" ")) ) );
    parameter_value = parameter_value.substr(parameter_value.find(" ")+1, parameter_value.length());
  }
  vector_to_deploy.push_back(stod( parameter_value ));
}

void Parameters::deploy_vector(vector<string> & vector_to_deploy, string parameter_value){
  vector_to_deploy.clear();
  while (parameter_value.find(" ") != string::npos){ 
    vector_to_deploy.push_back( parameter_value.substr(0, parameter_value.find(" ")) );
    parameter_value = parameter_value.substr(parameter_value.find(" ")+1, parameter_value.length());
  }
  vector_to_deploy.push_back(parameter_value );
}


void Parameters::deploy_matrix(vector<vector<double>> & matrix_to_deploy, string parameter_value){
  vector<double> vector_to_deploy;
  while (parameter_value.find(";") != string::npos){ 
    vector_to_deploy.clear();
    deploy_vector(vector_to_deploy, parameter_value.substr(0, parameter_value.find(";")));
    matrix_to_deploy.push_back(vector_to_deploy);
    parameter_value = parameter_value.substr(parameter_value.find(";")+2, parameter_value.length());
  }
  vector_to_deploy.clear();
  deploy_vector(vector_to_deploy, parameter_value);
  matrix_to_deploy.push_back(vector_to_deploy);
}
