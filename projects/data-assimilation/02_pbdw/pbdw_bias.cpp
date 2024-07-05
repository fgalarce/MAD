/*=============================================================================
  This file is part of the code MAD
  Copyright (C) 2017 - 2024,
    
     Felipe Galarce at WIAS

  MAD is free software; you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by the Free
  Software Foundation; either version 2.1 of the License, or (at your option)
  any later version.

  MAD is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with MDMA. If not, see http://www.gnu.org/licenses/.
  =============================================================================*/

#include <pbdw.hpp>

int main(int argc, char *argv[]){

  Parameters par = MADinitialize(argc, argv);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  PetscErrorCode code;

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);

  Geometry geo;
  geo.initialize(par, io);

  Boundary bd;
  bd.initialize(par, geo);

  InnerProduct ip;
  ip.initialize(par, geo, bd);

  LinearAlgebra la;
  la.initialize(par, ip);

  Calculus calculus;
  calculus.initialize(par, geo, bd);

  ROM rom;
  rom.initialize(par, geo, io, calculus);

  Observer observer;
  observer.initialize(par, geo, ip, io, calculus, la, bd);

//  observer.printRieszRepresenters();

  Mat G = buildCrossGramian(par, rom, observer, ip);
//  Mat G = buildCrossGramian(par, rom, observer, ip, par.dirResults() + "/G.bin");

  vector<double> sv = getSingularValues(G, par.nbModes());
  exportData(par.dirResults() + "/G_sv.txt", sv); 

  ofstream error_file(par.dirResults() + "/rec_error.txt");

  for (int t = par.start(); t < par.end(); t++){

    observer.buildMeasures(t);

    /* print observation on FEM basis in ensight */
    printObservations(par, observer, io, t);

    MADprint("Computing least-squares\n");
    Vec c = vec(par.nbModes());
    Mat N = mat(par.nbModes(), par.nbModes());
    MatMatMult(transpose(G), G, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &N);
    KSP ksp;
    configureKSP(ksp, par);
    KSPSetOperators(ksp, N, N);
    Vec rhs = vec(par.nbModes());
    MatMult(transpose(G), observer.measures(), rhs);
    KSPSolve(ksp, rhs, c);

    /* recover reconstruction from coordinates in ROM */
    vector<double> rom_coefficients = fingerprint(par, c, "./fingerprint_l2.txt");
//    vector<double> rom_coefficients = fingerprint(par, c);
    
    Vec u_tilde = vec(par.nbVariables()*io.nbVertices()); /* reconstruction */
    VecZeroEntries(u_tilde);
    for (int i = 0; i < par.nbModes(); i++){
      code = VecAXPY(u_tilde, rom_coefficients[i], rom.basis(i)); CHKERR(code);
    }
    for (int i = 0; i < par.nbVariables(); i++){
      io.writeState(calculus.split(u_tilde)[i], par.variableName()[i] + "_tilde", t);
    }

    /* \Pi_W u_tilde */
    Vec omega_tilde = zeros(par.nbVariables()*io.nbVertices());
    for (int i = 0; i < observer.nbMeasures(); i++){
      VecAXPY(omega_tilde, ip(observer.basis_l2(i), u_tilde), observer.basis_l2(i));
    }
    for (int i = 0; i < par.measureIt().size(); i++){
      io.writeState(calculus.split(omega_tilde)[i], par.variableName()[i] + "_omega_tilde", t);
    }

    /* l = 2*l_tilde - l_E(R(u)) */
    Vec theta_exp = observer.loadExpected("/Users/fgalarce/home/01_research/MAD/data/measures/MRI_2x2x2/IM_FFE_V60_0003_exp.txt");
    Vec new_coord(observer.measures());
    VecScale(new_coord, 2.0);
    VecAXPY(new_coord, -1.0, theta_exp);
   
    /* compute bias correction */ 
    Vec eta = vec(par.nbModes());
    MatMult(transpose(G), new_coord, eta);
    KSPSolve(ksp, eta, c);
   
    Vec u_star = vec(par.nbVariables()*io.nbVertices()); 
    VecZeroEntries(u_star);
    rom_coefficients = fingerprint(par, c, "./fingerprint_l2.txt");
//    rom_coefficients = fingerprint(par, c);
    for (int i = 0; i < par.nbModes(); i++){
      code = VecAXPY(u_star, rom_coefficients[i], rom.basis(i)); CHKERR(code);
    }
    for (int i = 0; i < par.nbVariables(); i++){
      io.writeState(calculus.split(u_star)[i], par.variableName()[i] + "_star", t);
    }

    /* compute error */
    vector<double> error_tilde = computeError(par, calculus, io, u_tilde, ip, t);
    vector<double> error_star = computeError(par, calculus, io, u_star, ip, t);
   
    for (int i = 0; i < par.nbVariables(); i++){ 
      MADprint("Relative error in var " + par.variableName()[i] + "0 = ", error_tilde[i]);
      MADprint("Relative error in var " + par.variableName()[i] + "* = ", error_star[i]);
      error_file << error_tilde[i] << " " << error_star[i] << " ";
    }
    error_file << endl;
  }

  MADfinalize(par);
}
