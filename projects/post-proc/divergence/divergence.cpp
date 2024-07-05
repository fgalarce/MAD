/*=============================================================================
  This file is part of the code M.D.M.A. (or MDMA)
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2020,
    
     Felipe Galarce at INRIA

  MDMA is free software; you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by the Free
  Software Foundation; either version 2.1 of the License, or (at your option)
  any later version.

  MDMA is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with MDMA. If not, see http://www.gnu.org/licenses/.
  =============================================================================*/

#include <ultra-4d-flow.hpp>

int main(int argc, char *argv[]){

  assert(argc == 2);

  string data_file = argv[1];
  
  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;
 
  Parameters par(data_file);
  par.print();

  /* Initialize MDMA */
  IO io;
  io.initialize(par);
    
  Geometry geo;
  geo.initialize(par, io);

  CFD cfd;
  cfd.initialize(par, geo, io);

  MasterTetrahedron tet;
  tet.initialize(par);

  InnerProduct ip;

  int nbVertices = cfd.io.nbVertices();
  int nbDofs = 3 * cfd.io.nbVertices(); 

  Tic tic;
  Mat mass = mat(nbVertices, nbVertices);
  int tetraId = 0; 
  for (vector<int> tetra : geo.tetrahedron()[0]){ /* loop on tetra */
    /* get finite element coordinates */
    vector<vector<double>> coordinates(4);
    coordinates[0] = geo.coordinates()[tetra[0]];
    coordinates[1] = geo.coordinates()[tetra[1]];
    coordinates[2] = geo.coordinates()[tetra[2]];
    coordinates[3] = geo.coordinates()[tetra[3]];

    tet.setCoordinates(coordinates);

    if (tetraId % (geo.tetrahedron()[0].size() / 10) == 0){
      cout << "  Computing elementary arrays for tetrahedron " << tetraId << " of " << geo.tetrahedron()[0].size() - 1 << endl;
      cout << "  Labels: " << tetra << endl; 
      cout << "  det(J) = " << tet.detJacobian() << endl; 
    }
    /* Assemble elementary matrices */
    for (int i = 0; i < 4; i++){
      for (int j = 0; j < 4; j++){
        /* mass */
        code = MatSetValue(mass, tetra[i], tetra[j], tet.mass(i,j), ADD_VALUES); CHKERRQ(code);
      }
    }
    tetraId++;
  }
  cout << " Time elapsed: " << tic.toc() << " sec. "<< endl;
  code = MatAssemblyBegin(mass, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(mass, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

  ofstream div_file(par.dirResults() + "/norm_div_u.txt");

  for (int iteration : range(par.start(), par.end())){

    cout << "\n - - - - - - - - - - - - - - - - - - -" << endl;
    cout << "Iteration: " << iteration << endl;
    cout << "Computing divergence." << endl;

    Vec u = vec(nbDofs);
    loadVec(u, par.dirSyntheticField() + "/" + par.patientName() + "." + wildcard(iteration) + ".vct");

    Vec divu = vec(nbVertices);
    divu = cfd.divergence(u);

    io.writeState(divu, "divu", iteration);

    double norm_divu = sqrt(ip(divu, mass, divu));
    cout << "L2-norm of div u = " << norm_divu << endl;

    div_file << norm_divu << endl;
  }

  div_file.close();
  SlepcFinalize();
}
