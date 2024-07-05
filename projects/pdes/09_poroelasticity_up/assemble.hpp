#include <mad.hpp>

void assembleSystem(Mat A, Mat  M_sol_2, Parameters par, Geometry geo, MasterElement fe){

  PetscErrorCode code;
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  int nbDofsPerNode = par.nbDofsPerNode()[0]; 
  int nbVertices = geo.nbVertices;
  int nbDofsSol = nbVertices * nbDofsPerNode;
  int nbDofsPress = nbVertices;
  int nbDofs = nbDofsSol + nbDofsPress;

  /* local values */
  int nbNodesPerElement = geo.dimension()+1; 
  double K_darcy =  par.permeability()/par.viscosity();

  if (world_rank == 0) cout << "PORO: permeability/viscosity = " << K_darcy << endl;

  double lambda = par.youngModulus()[1] * par.poissonRatio()[1] / ((1.0 + par.poissonRatio()[1])*(1.0 - 2.0*par.poissonRatio()[1]));
  double mu = par.youngModulus()[1]/(2.0*(1.0+par.poissonRatio()[1]));

  int m,n;
  code = MatGetOwnershipRange(A, &m, &n); CHKERR(code);
  if (world_rank == 0) cout << "PORO: Assembling discretization matrix." << endl;
  for (int partId = 0; partId < geo.elements().size(); partId++){

    if (world_rank == 0) cout << "  Assembling elementary matrices, part " << partId << endl;
    if (world_rank == 0) cout << "  mu = " << mu << endl;
    if (world_rank == 0) cout << "  lambda = " << lambda << endl;
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
              double Bij = fe.mixed(i,j,comp);
              code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp,  lambda*div_div_ij + mu*Kij, ADD_VALUES); CHKERR(code);
              code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp, nbDofsSol + simplex[j], Bij, ADD_VALUES); CHKERR(code);  /* -B */ 
            }
          }

          /* Darcy law for fluid  */
          if (nbDofsSol + simplex[i] >= m && nbDofsSol + simplex[i] < n){
            double Kij = fe.stiffness(i,j);
            code = MatSetValue(A, nbDofsSol + simplex[i], nbDofsSol + simplex[j], K_darcy * Kij, ADD_VALUES);
          }
          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
            if (nbDofsSol + simplex[i] >= m && nbDofsSol + simplex[i] < n){
              double Bij = fe.mixed(i,j,comp);
              code = MatSetValue(M_sol_2, nbDofsSol + simplex[i], nbDofsPerNode*simplex[j]+comp, Bij/par.timeStep(), ADD_VALUES); CHKERR(code); /* - B^T */
            }
          }
        }
      }
    }
  }
  code = MatAssemblyBegin(M_sol_2, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(M_sol_2, MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERR(code);

}
