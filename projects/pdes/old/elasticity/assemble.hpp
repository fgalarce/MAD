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

  double lambda = par.lambda();
  double mu = par.mu();

  if (world_rank == 0) cout << "Elasticity: mu = " << mu << " lambda = " << lambda << endl;

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
      for (int i = 0; i < nbNodesPerElement; i++){
        for (int j = 0; j < nbNodesPerElement; j++){
          /* - - - - - -
            Newton law for solid 
          - - - - - - - */
    
          /* \epsilon(u) : \epsilon(v)  */
          for (int comp = 0; comp < nbDofsPerNode; comp++){
            if (nbDofsPerNode*simplex[i]+comp >= m && nbDofsPerNode*simplex[i]+comp < n){
              double Kij = 1.0/2.0*fe.stiff(i,j); /* diagonal components */
              code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp,  2.0*mu*Kij, ADD_VALUES); CHKERR(code);
            }
          }
          for (int comp_v = 0; comp_v < nbDofsPerNode; comp_v++){
            if (nbDofsPerNode*simplex[i]+comp_v>= m && nbDofsPerNode*simplex[i]+comp_v < n){
              for (int comp_u = 0; comp_u < nbDofsPerNode; comp_u++){
                double Kdev = 1.0/2.0 * fe.dphi_dx(i, comp_u) * fe.dphi_dx(j, comp_v) * fe.detJacobian()*fe.weightsSum();
                code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp_v, nbDofsPerNode*simplex[j]+comp_u,  2.0*mu*Kdev, ADD_VALUES); CHKERR(code);
              }
            }
          }

          /* div(u) div(v) */
          for (int comp_v = 0; comp_v < nbDofsPerNode; comp_v++){
            /* deviatoric components */
            if (nbDofsPerNode*simplex[i]+comp_v>= m && nbDofsPerNode*simplex[i]+comp_v < n){
              for (int comp_u = 0; comp_u < nbDofsPerNode; comp_u++){
                double div_div_ij = fe.dphi_dx(i, comp_v) * fe.dphi_dx(j, comp_u) * fe.detJacobian()*fe.weightsSum();
                code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp_v, nbDofsPerNode*simplex[j]+comp_u,  lambda*div_div_ij, ADD_VALUES); CHKERR(code);
              }
            }
          }
        }
      }
    }
  }
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERR(code);

}
