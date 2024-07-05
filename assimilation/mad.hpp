/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2021,
    
     Felipe Galarce at INRIA / WIAS

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

#ifndef MAD
#define MAD

#include <iostream>
#include <fstream>
#include <vector>
#include <assert.h>
#include <cmath>

#include <slepc.h>
#include <STLvectorUtils.hpp>
#include <petscWrapper.hpp>
#include <tools.hpp>

#include <io.hpp>
#include <parameters.hpp>
#include <linearAlgebra.hpp>

/* fem */
#include <innerProduct.hpp>
#include <masterElement.hpp>
#include <masterElementBD.hpp>
#include <interpolate.hpp>
#include <boundaries.hpp>

#include <model.hpp>
#include <geometry.hpp>
//#include <cfd.hpp>
#include <measures.hpp>
#include <observer.hpp>
#include <rom.hpp>
#include <kalman.hpp>
#include <learning.hpp>
#include <calculus.hpp>

/* PDEs */
#include <laplacian.hpp>
#include <elasticity.hpp>
#include <elastodynamics.hpp>
//#include <ped_up.hpp>
//#include <pevd_up.hpp>
#include <navier_stokes.hpp>
#include <heat.hpp>
#include <wave.hpp>
#include <wave_complex.hpp>
#include <hammer.hpp>
//#include <natural_convection.hpp>
#include <poro_complex.hpp>
#include <ns.hpp>

/* Data assimilation */
#include <inverse_problems_tools.hpp>

using namespace std;

#define PI 3.141592 
#define EPSILON 0.001
#define GOOGOL 10000000000000

Parameters MADinitialize(int argc, char *argv[]){
  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  assert(argc == 2);
  /* Parse data file */
  string data_file = argv[1];
  Parameters par(data_file);
  par.print();
  return par;
}

void MADfinalize(Laplacian & laplacian, Parameters & par, Mat K, Vec rhs, Vec u){
  VecDestroy(&rhs);
  VecDestroy(&u);
  MatDestroy(&K);
  laplacian.finalize();
  par.finalize();
  SlepcFinalize();
}


void MADfinalize(Elastodynamics & elastodynamics, Mat A, Vec b, Vec u, Parameters & par){
  VecDestroy(&b);
  VecDestroy(&u);
  MatDestroy(&A);
  elastodynamics.finalize();
  par.finalize();
  SlepcFinalize();
}

void MADfinalize(Elasticity & elasticity, Mat A, Vec b, Vec u, Parameters & par){
  VecDestroy(&b);
  VecDestroy(&u);
  MatDestroy(&A);
  elasticity.finalize();
  par.finalize();
  SlepcFinalize();
}

void MADfinalize(Wave & wave, Mat A, Vec b, Vec u, Parameters & par){
  VecDestroy(&b);
  VecDestroy(&u);
  MatDestroy(&A);
  wave.finalize();
  par.finalize();
  SlepcFinalize();
}

void MADfinalize(WaveComplex & wave, Mat A, Vec b, Vec u, Parameters & par){
  VecDestroy(&b);
  VecDestroy(&u);
  MatDestroy(&A);
  wave.finalize();
  par.finalize();
  SlepcFinalize();
}

void MADfinalize(Hammer & hammer, Parameters & par){
  hammer.finalize();
  par.finalize_short();
  SlepcFinalize();
}

void MADfinalize(Heat & heat, Parameters & par, Mat K, Vec rhs, Vec u){
  VecDestroy(&rhs);
  VecDestroy(&u);
  MatDestroy(&K);
  heat.finalize();
  par.finalize();
  SlepcFinalize();
}

void MADfinalize(Parameters par){
  par.finalize();
  SlepcFinalize();
}

template <class generic>
void MADprint(string message, generic number) {
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank == 0) cout << message << number << endl;
}

void MADconsole(const char* cmd) {
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank == 0) system(cmd);
}

void MADprint(string message){
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if (world_rank == 0) cout << message;
}

#endif
