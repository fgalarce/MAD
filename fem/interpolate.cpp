/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2024,
    
     Felipe Galarce at INRIA / WIAS / PUCV

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

#include<interpolate.hpp>

void INT::initialize(Parameters parameters, const Geometry & geometry){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "INT: Initializing." << endl;
  par = parameters;
  geo = geometry;
  m_dimension = geo.dimension();
  assert(m_dimension == 2 or m_dimension == 3);
  m_nbDofs = m_dimension*geo.nbVertices;
  m_nbVertices = geo.nbVertices;
  m_verbose = par.verbose();

  fe.initialize(par, geo.dimension());
  feBD.initialize(par, geo.dimension());
}

void INT::initialize(Parameters parameters, const Geometry & geometry, const IO & inputOutput){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "INT: Initializing." << endl;
  par = parameters;
  geo = geometry;
  m_dimension = geo.dimension();
  assert(m_dimension == 2 or m_dimension == 3);
  m_nbDofs = m_dimension*geo.nbVertices;
  m_nbVertices = geo.nbVertices;
  m_verbose = par.verbose();
  io = inputOutput;
  fe.initialize(par, geo.dimension());
  feBD.initialize(par, geo.dimension());
}

Vec INT::interpolateP0_P1(vector<double> u){
  int nbDofsPerNode = 1;
  if (u.size() == geo.dimension()*geo.elements()[0].size()){
    nbDofsPerNode = geo.dimension();
  } 
  if (m_world_rank == 0) cout << "INT: Interpolating field to P1" << endl;
  int nbNodesPerElement = geo.dimension()+1;
  int nbElements = geo.elements()[0].size(); 
  Vec Iv = zeros(m_nbVertices*nbDofsPerNode); /* Iv = v / a point-wise */
  Vec v = zeros(m_nbVertices*nbDofsPerNode);
  Vec a = zeros(m_nbVertices*nbDofsPerNode);
  for (int feId = 0; feId < geo.elements()[0].size(); feId++){ /* loop on boundary elements */
    vector<int> simplex = geo.elements()[0][feId];
    /* get finite element coordinates */
    vector<vector<double>> coordinates(nbNodesPerElement);
    for (int i = 0; i < nbNodesPerElement; i++){
      coordinates[i] = geo.coordinates()[simplex[i]];
    }
    fe.setCoordinates(coordinates);
    fe.computeVolume();
    for (int i = 0; i < nbNodesPerElement; i++){
      for (int comp = 0; comp < nbDofsPerNode; comp++){
        vecSet(v, simplex[i] + comp*m_nbVertices, fe.volume() * u[feId + comp*nbElements]);
        vecSet(a, simplex[i] + comp*m_nbVertices, fe.volume());
      }
    }
  }

  VecAssemblyBegin(a);
  VecAssemblyEnd(a);
  VecAssemblyBegin(v);
  VecAssemblyEnd(v);

  VecPointwiseDivide(Iv, v, a);
  return Iv;
}


Vec INT::interpolateP0_P1_boundary(vector<double> u, string mode){
  int nbDofsPerNode = 1;
  int nbElementsBD = geo.nbElementsBoundary(); 
  if (u.size() == geo.dimension()*nbElementsBD){
    nbDofsPerNode = geo.dimension();
  } 
  if (m_world_rank == 0) cout << "INT: Interpolating boundary field to P1" << endl;
  int nbNodesPerElement = geo.dimension();

  MasterElementBD feBD;
  feBD.initialize(par, geo.dimension());
  int offset = 0;
  Vec Iv_full = zeros(m_nbVertices*nbDofsPerNode); /* Iv = v / a point-wise */
  for (int partId = 0; partId < geo.elementsBD().size(); partId++){ /* loop on parts */
    Vec Iv = zeros(m_nbVertices*nbDofsPerNode); /* Iv = v / a point-wise */
    Vec v = zeros(m_nbVertices*nbDofsPerNode);
    Vec a = zeros(m_nbVertices*nbDofsPerNode);
    for (int feId = 0; feId < geo.elementsBD()[partId].size(); feId++){ /* loop on boundary elements */
      vector<int> simplex = geo.elementsBD()[partId][feId];
      /* get finite element coordinates */
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int i = 0; i < nbNodesPerElement; i++){
        coordinates[i] = geo.coordinates()[simplex[i]];
      }
      feBD.setCoordinates(coordinates);

      feBD.computeVolume();
      for (int i = 0; i < nbNodesPerElement; i++){
        for (int comp = 0; comp < nbDofsPerNode; comp++){
          if (mode == "splitted"){
            vecSet(v, simplex[i] + comp*m_nbVertices, feBD.volume() * u[offset + feId + comp*nbElementsBD]);
            vecSet(a, simplex[i] + comp*m_nbVertices, feBD.volume());
          } else if (mode == "contiguous"){
            vecSet(v, m_dimension*simplex[i] + comp, feBD.volume() * u[offset + m_dimension*feId + comp]);
            vecSet(a, m_dimension*simplex[i] + comp, feBD.volume());
          }
        }
      }
    }
    if (mode == "splitted"){
      offset = offset + geo.elementsBD()[partId].size();
    } else if (mode == "contiguous"){
      offset = offset + geo.elementsBD()[partId].size()*m_dimension;
    }
    VecAssemblyBegin(v);
    VecAssemblyEnd(v);
    VecAssemblyBegin(a);
    VecAssemblyEnd(a);

    /* fill zero elements */
    vector<double> aSTL = stl(a);
    for (int i = 0; i < nbDofsPerNode*m_nbVertices; i++){
      if (aSTL[i] == 0.0){
        vecSetInsert(a, i, 1.0);
      }
    }
    VecAssemblyBegin(a);
    VecAssemblyEnd(a);

    VecPointwiseDivide(Iv, v, a);
    vector<double> Iv_stl = stl(Iv);
    /* This avoids problems in the corners when too sharp, but it does depend on the label order, eventually it will lead to bugs */
    for (int i = 0; i <  Iv_stl.size(); i++){
      if (Iv_stl[i] != 0.0){
        vecSetInsert(Iv_full, i, Iv_stl[i]);
      }
    }
  }

  return Iv_full;
}

Vec INT::interpolate_field(string filePath, IO & io_i, int varLabel){
  int nbVertices_i = io_i.nbVertices();
  int nbDofs_i = par.nbDofsPerNode()[varLabel]*nbVertices_i;
  Vec u_i = vec(nbDofs_i);
  io_i.loadVector(u_i, filePath);
  return interpolate_field(u_i, io_i);
}

Vec INT::interpolate_field(Vec u_i, IO & io_i, int varLabel){

  if (m_world_rank == 0 ) cout << "Interpolation: P1 barycentric interpolation." << endl;

  Geometry geo_i;
  geo_i.initialize(io_i);     

  int nbVertices_i = io_i.nbVertices();
  int nbVertices_j = io.nbVertices();
  int nbDofs_i = par.nbDofsPerNode()[varLabel]*nbVertices_i;
  int nbDofs_j = par.nbDofsPerNode()[varLabel]*nbVertices_j; 
  int nbDofsPerNode = par.nbDofsPerNode()[varLabel];
  int nbNodesPerElement = geo_i.dimension() + 1;

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        Compute I u (\Omega_i) 
        \Omega_j : the fields will be interpolated here 
      - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /* get stl version of the vector (no way to paralelize this) */
  vector<double> u(nbDofs_i, 0.0);
  Vec u_seq = getSequential(u_i);
  u = stl(u_seq);

  int low, high;
  Vec Iu = zeros(nbDofs_j);
  VecGetOwnershipRange(Iu, &low, &high);

  bool * found = new bool[nbDofsPerNode*geo.coordinates().size()];
  for (int k=0; k<geo.coordinates().size(); k++){
    for (int comp=0; comp<nbDofsPerNode; comp++){
      found[nbDofsPerNode*k+comp] = false;
    }
  }

  MasterElement fe;
  fe.initialize(par, geo.dimension());

  for (int k = 0; k < geo.coordinates().size(); k++){
      
    if (k % (geo.coordinates().size() / par.verbose()) == 0){
      if (m_world_rank == 0) cout << "  Vertex " << k << " / " <<  geo.coordinates().size() << endl;
    }

    for (int comp = 0; comp < nbDofsPerNode; comp++){
      if( nbDofsPerNode*k + comp >= low && nbDofsPerNode*k + comp < high){
        for (vector<int> simplex : geo_i.elements()[0]){
          /* Set finite element */
          vector<vector<double>> coord(geo_i.dimension()+1);
          for (int idNodeCoord = 0; idNodeCoord < geo_i.dimension() + 1; idNodeCoord++){
            coord[idNodeCoord] = geo_i.coordinates()[simplex[idNodeCoord]];
          }
          fe.setCoordinates(coord);
          fe.computeSize();
          double epsilon = fe.size()*1e-9;
          int flagBC = 0;
          vector<double> bc = fe.barycentric_coor(geo.coordinates()[k]);
          for (int checkBC = 0; checkBC < geo.dimension() + 1; checkBC++){
            if (bc[checkBC] >= -epsilon){
              flagBC++;
            }
          }
          if (flagBC = geo.dimension()+1){
            double nodalValue = 0.0;
            for (int idNode = 0; idNode < nbNodesPerElement; idNode++){
              nodalValue = nodalValue + u[nbDofsPerNode*simplex[idNode] + comp] * bc[idNode];
            }
            for (int idNode = 0; idNode < nbNodesPerElement; idNode++){
              code = VecSetValue(Iu, nbDofsPerNode*k + comp, nodalValue, INSERT_VALUES); CHKERR(code);
            }
            found[nbDofsPerNode*k + comp] = true;
            code = MPI_Bcast(&(found[nbDofsPerNode*k + comp]), 1, MPI_CXX_BOOL, m_world_rank, MPI_COMM_WORLD); CHKERR(code);
            break;
          }
        }
      }
    }
  }
  code = MPI_Barrier(MPI_COMM_WORLD); CHKERR(code);

  if (m_world_rank == 0 ) cout << "Interpolation: correcting via NN." << endl;
  for (int k = 0; k < geo.coordinates().size(); k++){
      
    if (k % (geo.coordinates().size() / par.verbose()) == 0){
      if (m_world_rank == 0) cout << "  Vertex " << k << " / " <<  geo.coordinates().size() << endl;
    }
    for (int comp = 0; comp < nbDofsPerNode; comp++){
      if (!found[nbDofsPerNode*k+comp]){
        if( nbDofsPerNode*k + comp >= low && nbDofsPerNode*k + comp < high ){
          int nn = geo_i.nearest_neighbourd(geo.coordinates()[k]);
          code = VecSetValue(Iu, nbDofsPerNode*k + comp, u[nbDofsPerNode*nn + comp], INSERT_VALUES); CHKERR(code);
        }
      }
    }
  }
  code = VecAssemblyBegin(Iu); CHKERR(code);
  code = VecAssemblyEnd(Iu); CHKERR(code);

  return Iu;
}

vector<Vec> INT::interpolate_field(vector<Vec> u_i, IO & io_i, int varLabel){

  Geometry geo_i;
  geo_i.initialize(io_i);     

  int nbVertices_i = io_i.nbVertices();
  int nbVertices_j = io.nbVertices();
  int nbDofs_i = par.nbDofsPerNode()[varLabel]*nbVertices_i;
  int nbDofs_j = par.nbDofsPerNode()[varLabel]*nbVertices_j; 
  int nbDofsPerNode = par.nbDofsPerNode()[varLabel];
  int nbNodesPerElement = geo_i.dimension() + 1;

  /*  - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        Compute I u (\Omega_i) 
        \Omega_j : the fields will be interpolated here 
      - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /* get stl version of the vector (no way to paralelize this) */
  int low, high;
  vector<Vec> Iu(u_i.size());
  vector<vector<double>> u(u_i.size());
  for (int i = 0; i < u_i.size(); i++){
    Vec u_seq = getSequential(u_i[i]);
    u[i].resize(nbDofs_i);
    u[i] = stl(u_seq);
    vec(Iu[i], nbDofs_j);
    VecZeroEntries(Iu[i]);
  }

  code = VecGetOwnershipRange(Iu[0], &low, &high); CHKERR(code);

  bool * found = new bool[nbDofsPerNode*geo.coordinates().size()];
  for (int k=0; k<geo.coordinates().size(); k++){
    for (int comp=0; comp<nbDofsPerNode; comp++){
      found[nbDofsPerNode*k+comp] = false;
    }
  }

  vector<vector<int>> elements;
  for (int idRegion = 0; idRegion < geo_i.elements().size(); idRegion++){
    for (int idEl = 0; idEl < geo_i.elements()[idRegion].size(); idEl++){
      elements.push_back(geo_i.elements()[idRegion][idEl]);
    }
  }

  MasterElement fe;
  fe.initialize(par, geo.dimension());

  if (m_world_rank == 0 ) cout << "Interpolation: running P1 approximation." << endl;

  for (int k = 0; k < geo.coordinates().size(); k++){
    if (k % (geo.coordinates().size() / par.verbose()) == 0){
      if (m_world_rank == 0) cout << "  Vertex " << k << " / " <<  geo.coordinates().size() << endl;
    }

    for (int comp = 0; comp < nbDofsPerNode; comp++){
      if( nbDofsPerNode*k + comp >= low && nbDofsPerNode*k + comp < high){
        for (vector<int> simplex : elements){
          /* Set finite element */
          vector<vector<double>> coord(geo_i.dimension()+1);
          for (int idNodeCoord = 0; idNodeCoord < geo_i.dimension() + 1; idNodeCoord++){
            coord[idNodeCoord] = geo_i.coordinates()[simplex[idNodeCoord]];
          }
          fe.setCoordinates(coord);
          vector<double> bc = fe.barycentric_coor(geo.coordinates()[k]);

          if ( bc[0] > 0 && bc[1] > 0 && bc[2] > 0 ){
            for (int i = 0; i < u_i.size(); i++){
              double nodalValue = 0.0;
              for (int idNode = 0; idNode < nbNodesPerElement; idNode++){
                nodalValue = nodalValue + u[i][nbDofsPerNode*simplex[idNode] + comp] * bc[idNode];
              }
              for (int idNode = 0; idNode < nbNodesPerElement; idNode++){
                code = VecSetValue(Iu[i], nbDofsPerNode*k + comp, nodalValue, INSERT_VALUES); CHKERR(code);
              }
            }
            found[nbDofsPerNode*k + comp] = true;
            code = MPI_Bcast(&(found[nbDofsPerNode*k + comp]), 1, MPI_CXX_BOOL, m_world_rank, MPI_COMM_WORLD); CHKERR(code);
            break;
          }
        }
      }
    }
  }
  code = MPI_Barrier(MPI_COMM_WORLD); CHKERR(code);

  if (m_world_rank == 0 ) cout << "Interpolation: correcting via NN." << endl;
  for (int k = 0; k < geo.coordinates().size(); k++){
      
    if (k % (geo.coordinates().size() / par.verbose()) == 0){
      if (m_world_rank == 0) cout << "  Vertex " << k << " / " <<  geo.coordinates().size() << endl;
    }
  
    for (int comp = 0; comp < nbDofsPerNode; comp++){
      if (!found[nbDofsPerNode*k+comp]){
        if( nbDofsPerNode*k + comp >= low && nbDofsPerNode*k + comp < high ){
          int nn = geo_i.nearest_neighbourd(geo.coordinates()[k]);
          for (int i = 0; i < u_i.size(); i++){
            code = VecSetValue(Iu[i], nbDofsPerNode*k + comp, u[i][nbDofsPerNode*nn + comp], INSERT_VALUES); CHKERR(code);
          }
        }
      }
    }
  }

  for (int i = 0; i < u_i.size(); i++){
    code = VecAssemblyBegin(Iu[i]); CHKERR(code);
    code = VecAssemblyEnd(Iu[i]); CHKERR(code);
  }

  return Iu;
}
