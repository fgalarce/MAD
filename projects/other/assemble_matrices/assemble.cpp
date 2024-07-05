#include<ultra-4d-flow.hpp>
#include<cfd.hpp>

int main(int argc, char *argv[]){

  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  /* Parse data file */
  assert(argc == 2);
  string data_file = argv[1];
  Parameters par(data_file);
  par.print();

//  CFD cfd;
//  cfd.initialize(par);

  /* Initialize MDMA objects */
  IO io;
  io.initialize(par);

  Geometry geo;
  geo.initialize(par, io);
  
  int nbDofs = io.nbVertices();
  int nbVertices = io.nbVertices();

  int nbDofsPerNode = 1;
  if (par.problemType() == "vector"){
    nbDofs = nbDofs * 3;
    nbDofsPerNode = nbDofsPerNode * 3;
  }

  /* Initialize finite element */
  MasterTetrahedron tetrahedron;
  tetrahedron.initialize(par);

  Mat mass = mat(nbDofs, nbDofs);
  Mat stiff = mat(nbDofs, nbDofs);

  Tic tic;
  cout << "Assembling mass and stiffness matrix\n";
  for (int partId = 0; partId < geo.tetrahedron().size(); partId++){ /* loop on parts */

    for (vector<int> tetra : geo.tetrahedron()[partId]){ /* loop on tetra */
      /* get finite element coordinates */
      vector<vector<double>> coordinates(4);
      coordinates[0] = geo.coordinates()[tetra[0]];
      coordinates[1] = geo.coordinates()[tetra[1]];
      coordinates[2] = geo.coordinates()[tetra[2]];
      coordinates[3] = geo.coordinates()[tetra[3]];

      tetrahedron.setCoordinates(coordinates);
        
      /* set values from local to global */
      for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
          for (int comp = 0; comp < nbDofsPerNode; comp++){
            code = MatSetValue(mass, nbDofsPerNode * tetra[i] + comp, nbDofsPerNode * tetra[j] + comp, tetrahedron.mass(i,j), ADD_VALUES); CHKERRQ(code);
            code = MatSetValue(stiff, nbDofsPerNode *  tetra[i] + comp, nbDofsPerNode * tetra[j] + comp, tetrahedron.stiffness(i,j), ADD_VALUES); CHKERRQ(code);
          }
        }
      }
    }
  }
  cout << "Time elapsed: " << tic.toc() / 60.0 << " minutes." << endl;

  code = MatAssemblyBegin(mass, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(mass, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyBegin(stiff, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(stiff, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

  cout << "Norm mass: " << norm(mass) << endl;
  cout << "Norm stiff: " << norm(stiff) << endl;

  mkdir(par.dirResults());
  saveMat(mass, par.dirResults() + "/mass.bin");
  saveMat(stiff, par.dirResults() + "/stiffness.bin");

}
