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

  Calculus calculus;
  calculus.initialize(par, geo, bd);

  ROM rom;
  rom.initialize(par, geo, io, calculus);

  /* Computing number of snaps */
  int jump = 1;
  int nbSnapshots = 0;
  bool static_problem = false;
  if (static_problem){
    nbSnapshots = par.nbSimulations();
  } else {
    for (int simId = 0; simId < par.nbSimulations(); simId++){
      if (par.snaps_per_sim() == -1){
        nbSnapshots = nbSnapshots + floor(stoi(getParameter("nbIterations", par.maniFolder() + "sim" + wildcard(simId) + "/parameters"))/jump) - par.start();
      } else {
        nbSnapshots = nbSnapshots + par.snaps_per_sim();
      }
    }
  }

  /* Snapshots */
  int nbDofs = 0;
  for (int i = 0; i < par.nbVariables(); i++){
    nbDofs = nbDofs + io.nbVertices()*par.nbDofsPerNode()[i];
  }
  vector<Vec> V(nbSnapshots);
  int m,n;
  vec(V[0], nbDofs);
  code = VecGetOwnershipRange(V[0], &m, &n); CHKERR(code);
  int nbDofsLocal = n - m;

  MADprint("POD: reading data. # of snapshots: ", nbSnapshots); 
  int iterSnap = 0;
  for (int simId = 0; simId < par.nbSimulations(); simId++){
    int nbSnaps;
    if (static_problem){
      nbSnaps = 1;
    } else {
      if (par.snaps_per_sim() == -1){
        nbSnaps = floor(stoi(getParameter("nbIterations", par.maniFolder() + "sim" + wildcard(simId) + "/parameters")));
      } else {
        nbSnaps = par.snaps_per_sim() + par.start();
      }
    }
    for (int i = par.start(); i < nbSnaps; i = i + jump){
      vec(V[iterSnap], nbDofs);
      vector<string> to_load(par.nbVariables());
      for (int var = 0; var < par.nbVariables(); var++){
        string end_str;
        if (par.nbDofsPerNode()[var] == 1){
          to_load[var] = par.maniFolder() + "sim" + wildcard(simId) + "/" + par.variableName()[var] + "." + wildcard(i) + ".scl"; 
        } else {
          to_load[var] = par.maniFolder() + "sim" + wildcard(simId) + "/" + par.variableName()[var] + "." + wildcard(i) + ".vct"; 
        }
      }
      io.loadVector(V[iterSnap], to_load);
      code = VecAssemblyBegin(V[iterSnap]); CHKERRQ(code);
      code = VecAssemblyEnd(V[iterSnap]); CHKERRQ(code);
      double snapshot[nbDofsLocal];
      code = VecGetValues(V[iterSnap], nbDofsLocal, &range(m, n)[0], snapshot); CHKERRQ(code); 
      iterSnap = iterSnap + 1;
    }
  }

  MADprint("Computing fingerprint");
  vector<vector<double>> fingerprint(par.nbModes());
  for (int i = 0; i < par.nbModes(); i++){
    fingerprint[i].resize(2, 0.0);
    for (int j = 0; j < V.size(); j++){
      double fp = ip(V[j], rom.basis(i));
      if (fp < fingerprint[i][0]){
        fingerprint[i][0] = fp;
      }
      if (fp > fingerprint[i][1]){
        fingerprint[i][1] = fp;
      } 
    }
  }
  exportData(par.dirResults() + "/fingerprint.txt", fingerprint);

  MADfinalize(par);
}
