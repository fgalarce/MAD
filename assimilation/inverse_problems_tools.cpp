/*=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2023,
    
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

#include<inverse_problems_tools.hpp>

void InverseProblemsTools::initialize(Parameters parameters, const Geometry & geometry){
  MPI_Comm_rank(MPI_COMM_WORLD, &m_world_rank);
  if (m_world_rank == 0) cout << "InverseProblemsTools: Initializing" << endl;
  par = parameters;
  geo = geometry;
  fe.initialize(par, geo.dimension());
  m_verbose = par.verbose();
  nbDofs = 0;
  assert(par.nbDofsPerNode().size() == 1);

  nbDofsPerNode = par.nbDofsPerNode()[0]; 
  nbVertices = geo.nbVertices;
  nbDofs = nbVertices * nbDofsPerNode;

}

void InverseProblemsTools::finalize(){
  MatDestroy(&M);
  MatDestroy(&K);
}

void InverseProblemsTools::assembleMassAndStiffness(){

  mat(M, nbDofs, nbDofs); 
  mat(K, nbDofs, nbDofs); 

  /* local values */
  int nbNodesPerElement = geo.dimension()+1; 

  int m,n;
  code = MatGetOwnershipRange(K, &m, &n); CHKERR(code);

  for (int partId = 0; partId < geo.elements().size(); partId++){

    if (m_world_rank == 0) cout << "InverseProblemsTools: assembling mass and stiffness matrix. Part " << partId << endl;

    MasterElement fe; 
    fe.initialize(par, geo.dimension());

    for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on finiteElement */
      
      /* get finite element coordinates */
      vector<int> simplex = geo.elements()[partId][feId];
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int nodeId = 0; nodeId < nbNodesPerElement; nodeId++){
        coordinates[nodeId] = geo.coordinates()[simplex[nodeId]];
      }
      fe.setCoordinates(coordinates);

      if (feId % (geo.elements()[partId].size() / par.verbose()) == 0){
        if (m_world_rank == 0) cout << "    Elementary matrix for element: " << feId << "/" << geo.elements()[partId].size() - 1 << endl;
        if (m_world_rank == 0) cout << "    detJac: " << fe.detJacobian() << endl;
      }
      /* Assemble elementary matrices */
      for (int i = 0; i < nbNodesPerElement; i++){
        for (int j = 0; j < nbNodesPerElement; j++){

          /* mass */
          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
            if (nbDofsPerNode*simplex[i]+comp >= m && nbDofsPerNode*simplex[i]+comp < n){
              double Mij = fe.mass(i,j);
              code = MatSetValue(M, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp, Mij, ADD_VALUES); CHKERR(code);
            }
          }

          /* stifness  */
          for (int comp = 0; comp < nbDofsPerNode; comp++){
            if (nbDofsPerNode*simplex[i]+comp >= m && nbDofsPerNode*simplex[i]+comp < n){
              double Kij = fe.stiff(i,j); /* diagonal components */
              code = MatSetValue(K, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp, Kij, ADD_VALUES); CHKERR(code);
            }
          }
        }
      }
    }
  }

  code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY); CHKERR(code);
}
