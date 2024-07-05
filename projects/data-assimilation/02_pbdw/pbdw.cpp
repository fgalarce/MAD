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

#include <mad.hpp>

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

  for (int t = par.start(); t < par.end(); t++){

    /* build coordinates from real data */
    observer.buildMeasures(t);

    /* print synthetic observation in ensight */
    vector<Vec> obs = observer.observations(t);
    for (int i = 0; i < par.measureIt().size(); i++){
      io.writeState(obs[i], "measures_" + par.variableName()[i], t);
    }

    MADprint("Computing cross-Gramian\n");
    Mat G = mat(observer.nbMeasures(), par.nbModes());
    for (int i = 0; i < observer.nbMeasures(); i++){
      MADprint(to_string(float(i)/float(observer.nbMeasures())*100) + " %   \r");
      for (int j = 0; j < par.nbModes(); j++){
        matSetInsert(G, i, j, ip(observer.basis(i), rom.basis(j)));
      }
    }
    MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY);

    vector<double> svG = getSingularValues(G, par.nbModes());
    exportData("./svG.txt", svG);

    MADprint("Computing least-squares\n");
    Vec u = vec(par.nbVariables()*io.nbVertices()); /* reconstruction */
    Vec c = vec(par.nbModes());
    Mat N = mat(par.nbModes(), par.nbModes());
    MatMatMult(transpose(G), G, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &N);
    Vec rhs = vec(par.nbModes());
    MatMult(transpose(G), observer.measures(), rhs);
    KSP ksp;
    configureKSP(ksp, par);
    KSPSetOperators(ksp, N, N);
    KSPSolve(ksp, rhs, c);

    /* recover reconstruction from coordinates in ROM */
    VecZeroEntries(u);
    vector<double> c_stl = stl(c);
    for (int i = 0; i < par.nbModes(); i++){
      code = VecAXPY(u, c_stl[i], rom.basis(i)); CHKERR(code);
    }
    for (int i = 0; i < par.nbVariables(); i++){
      io.writeState(calculus.split(u)[i], par.variableName()[i] + "_star", t);
    }
  }

  MADfinalize(par);
}
