/*=============================================================================
  This file is part of the code MAD
  Copyright (C) 2017-2024,
    
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

#include<rom.hpp>

void ROM::initialize(Parameters parameters, Geometry & geometry, IO & inputOutput, Calculus & calc){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);

  par = parameters;
  geo = geometry;
  calculus = calc;
  io = inputOutput;

  m_nbDofs = 0;
  for (int i = 0; i < par.nbVariables(); i++){
    m_nbDofs = m_nbDofs + geo.nbVertices*par.nbDofsPerNode()[i];
  }

  /* load Basis */
  m_basis.resize(par.nbModes());
  for (int i = 0; i < par.nbModes(); i++){
    vec(m_basis[i], m_nbDofs);
    vector<string> to_load(par.nbVariables());
    for (int j = 0; j < par.nbVariables(); j++){
      to_load[j] = par.dirModel() + "/" + par.variableName()[j] + "_mode." + wildcard(i) + ".scl";
    }
    io.loadVector(m_basis[i], to_load);
  }
}
