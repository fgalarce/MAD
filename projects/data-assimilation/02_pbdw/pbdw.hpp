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

Mat buildCrossGramian(Parameters par, ROM & rom, Observer & observer, InnerProduct & ip, string dirCG = " ");
Mat assembleCrossGramian(Parameters par, ROM & rom, Observer & observer, InnerProduct & ip);
void printObservations(Parameters par, Observer & observer, IO & io, int t);
vector<double> fingerprint(Parameters par, Vec c, string fileFingerprint = " ");
double computeError(Parameters par, IO & io, Vec & u, InnerProduct & ip);

Mat buildCrossGramian(Parameters par, ROM & rom, Observer & observer, InnerProduct & ip, string dirCG){
  Mat G = mat(observer.nbMeasures(), par.nbModes());
  if (dirCG != " "){
    loadMat(G, dirCG);
  } else {
    G = assembleCrossGramian(par, rom, observer, ip);
  }
  return G;
}

Mat assembleCrossGramian(Parameters par, ROM & rom, Observer & observer, InnerProduct & ip){
  MADprint("Computing cross-Gramian matrix\n");
  Mat G = mat(observer.nbMeasures(), par.nbModes());
  for (int i = 0; i < observer.nbMeasures(); i++){
    MADprint(to_string(float(i)/float(observer.nbMeasures())*100) + " %   \r");
    for (int j = 0; j < par.nbModes(); j++){
      matSetInsert(G, i, j, ip(observer.basis(i), rom.basis(j)));
    }
  }
  MatAssemblyBegin(G, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(G, MAT_FINAL_ASSEMBLY);
  saveMat(G, par.dirResults() + "/G.bin");
  return G;
}

void printObservations(Parameters par, Observer & observer, IO & io, int t){
  vector<Vec> obs = observer.observations(t);
  for (int i = 0; i < par.measureIt().size(); i++){
    io.writeState(obs[i], "measures_" + par.variableName()[i], t);
  }
}

vector<double> fingerprint(Parameters par, Vec c, string fileFingerprint){
  vector<double> c_stl = stl(c);
  if (fileFingerprint != " "){
    vector<vector<double>> fin = importdata(fileFingerprint);
    for (int i = 0; i < par.nbModes(); i++){
      if (c_stl[i] < fin[i][0]){
        c_stl[i] = fin[i][0];
      }
      if (c_stl[i] > fin[i][1]){
        c_stl[i] = fin[i][1];
      }
    }
  }   
  return c_stl;
}

vector<double> computeError(Parameters par, Calculus & calculus, IO & io, Vec & u, InnerProduct & ip, int iteration){
  Vec u_gt = vec(io.nbVertices()*par.nbVariables());
  vector<string> to_load(par.nbVariables());
  for (int j = 0; j < par.nbVariables(); j++){
    to_load[j] = par.dirSyntheticField() + "/" + par.variableName()[j] + "." + wildcard(iteration) + ".scl";
  } 
  io.loadVector(u_gt, to_load);

  vector<double> error_rel(par.nbVariables());
  for (int i = 0; i < par.nbVariables(); i++){
    Vec error;
    zeros(error, io.nbVertices());
    VecAXPY(error, 1.0, calculus.split(u_gt)[i]);
    VecAXPY(error, -1.0, calculus.split(u)[i]);
    error_rel[i] = sqrt(dot(error, error)/dot(calculus.split(u_gt)[i], calculus.split(u_gt)[i]));
  }
  return error_rel;
//  return sqrt(ip(error, error)/ip(u_gt, u_gt));
}


double computeErrorIP(Parameters par, IO & io, Vec & u, InnerProduct & ip, int iteration){
  Vec u_gt = vec(io.nbVertices()*par.nbVariables());
  vector<string> to_load(par.nbVariables());
  for (int j = 0; j < par.nbVariables(); j++){
    to_load[j] = par.dirSyntheticField() + "/" + par.variableName()[j] + "." + wildcard(iteration) + ".scl";
  } 
  io.loadVector(u_gt, to_load);

  double error_rel(par.nbVariables());
  Vec error = vec(par.nbVariables()*io.nbVertices());
  VecZeroEntries(error);
  VecAXPY(error, 1.0, u_gt);
  VecAXPY(error, -1.0, u);
  error_rel = ip(error, error)/ip(u_gt, u_gt);
  return error_rel;
//  return sqrt(ip(error, error)/ip(u_gt, u_gt));
}
