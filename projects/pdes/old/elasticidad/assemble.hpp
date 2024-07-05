#include <mad.hpp>

void assembleSystem(Mat A, Parameters par, Geometry geo, MasterElement fe){

  PetscErrorCode code;
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  int nbDofsPerNode = par.nbDofsPerNode()[0]; 
  int nbVertices = geo.nbVertices;
  int nbDofs = nbVertices*nbDofsPerNode;

  /* local values */
  int nbNodesPerElement = geo.dimension()+1; 

  double lambda = par.youngModulus()[0] * par.poissonRatio()[0] / ((1.0 + par.poissonRatio()[1])*(1.0 - 2.0*par.poissonRatio()[0]));
  double mu = par.youngModulus()[0]/(1.0+par.poissonRatio()[0]);

  int m,n;
  code = MatGetOwnershipRange(A, &m, &n); CHKERR(code);
  if (world_rank == 0) cout << "Elasticity: Assembling discretization matrix." << endl;
  for (int partId = 0; partId < geo.elements().size(); partId++){

    if (world_rank == 0) cout << "  Assembling elementary matrices, part " << partId << endl;
    for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on elements */

      vector<int> simplex = geo.elements()[partId][feId];
      /* get finite element coordinates */
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int i = 0; i < nbNodesPerElement; i++){
        coordinates[i] = geo.coordinates()[simplex[i]];
      }

      fe.setCoordinates(coordinates);

      if (feId % (geo.elements()[partId].size() / par.verbose()) == 0){
        if (world_rank == 0) cout << "    Elementary matrix for elemnt: " << feId << "/" << geo.elements()[partId].size() - 1 << endl;
      }

      /* Assemble elementary matrices */
      for (int j = 0; j < nbNodesPerElement; j++){
          
        for (int i = 0; i < nbNodesPerElement; i++){
          /* Newton law for solid */
          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
            if (nbDofsPerNode*simplex[i]+comp >= m && nbDofsPerNode*simplex[i]+comp < n){
              double Kij = fe.stiffness_symmetric(i,j);
              double div_div_ij = fe.div_div(i,j,comp);
              code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp,  lambda*div_div_ij + mu*Kij, ADD_VALUES); CHKERR(code);
            }
          }
        }
      }
    }
  }
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERR(code);

}
