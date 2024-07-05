#include <mad.hpp>

void assembleSystem(Mat M, Mat B, Parameters par, Geometry geo, MasterElement fe){
  
  PetscErrorCode code;
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  int nbDofs = (par.nbDofsPerNode()[0] + par.nbDofsPerNode()[1])*geo.nbVertices;
  int nbDofsPress = geo.nbVertices;
  int nbDofsVel = par.nbDofsPerNode()[0]*geo.nbVertices;
  int nbVertices = geo.nbVertices;

  /* local values */
  int nbDofsPerNode = par.nbDofsPerNode()[0]; 
  int nbNodesPerElement = geo.dimension()+1; 

  int m,n;
  code = MatGetOwnershipRange(M, &m, &n); CHKERR(code);
  if (world_rank == 0) cout << "NS: Assembling discretization matrix." << endl;
  for (int partId = 0; partId < geo.elements().size(); partId++){
    for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on elements */

      vector<int> simplex = geo.elements()[partId][feId];
      /* get finite element coordinates */
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int i = 0; i < nbNodesPerElement; i++){
        coordinates[i] = geo.coordinates()[simplex[i]];
      }

      fe.setCoordinates(coordinates);
      fe.computeSize();

      if (feId % (geo.elements()[partId].size() / 5) == 0){
        if (world_rank == 0) cout << "  Elementary matrix for tetra: " << feId << "/" << geo.elements().size() - 1 << endl;
      }
      /* Assemble elementary matrices */
      for (int j = 0; j < nbNodesPerElement; j++){
        for (int i = 0; i < nbNodesPerElement; i++){
          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
            if (par.nbDofsPerNode()[0]*simplex[i]+comp >= m && par.nbDofsPerNode()[0]*simplex[i]+comp < n){
              double Kij = fe.stiffness(i,j);
              double Mij = par.density()/par.timeStep()*fe.mass(i,j);
              /* stiffness */
              code = MatSetValue(M, par.nbDofsPerNode()[0]*simplex[i]+comp, par.nbDofsPerNode()[0]*simplex[j]+comp, Mij, ADD_VALUES); CHKERR(code);
              code = MatSetValue(B, par.nbDofsPerNode()[0]*simplex[i]+comp, par.nbDofsPerNode()[0]*simplex[j]+comp, par.viscosity() * Kij, ADD_VALUES); CHKERR(code);
              /* div_phi_i phi_j */
              double Bij = fe.mixed(i,j,comp);
              code = MatSetValue(B, par.nbDofsPerNode()[0]*simplex[i]+comp, nbDofsVel + simplex[j], -1.0*Bij, ADD_VALUES); CHKERR(code); /* -B*/ 
            }
            if (nbDofsVel + simplex[j] >= m && nbDofsVel + simplex[j] < n){
              double Bij = fe.mixed(i,j,comp);
              code = MatSetValue(B, nbDofsVel + simplex[j], par.nbDofsPerNode()[0]*simplex[i]+comp, 1.0*Bij, ADD_VALUES); CHKERR(code); /* B^T */ 
            }
          }
          if (nbDofsVel + simplex[i] >= m && nbDofsVel + simplex[i] < n){
            double Kij = fe.stiffness(i,j);
            /* Brezzi-Pitkaranta stabilization for P1-P1 elements */
            code = MatSetValue(B, nbDofsVel + simplex[i], nbDofsVel + simplex[j], fe.size()*fe.size() * Kij, ADD_VALUES); CHKERR(code);
          }
        }
      }
    }
  }
  code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY); CHKERR(code);

}

void assembleConvection(Mat C, Vec phi, Parameters par, Geometry geo, MasterElement fe){

    PetscErrorCode code;
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int nbDofs = (par.nbDofsPerNode()[0] + par.nbDofsPerNode()[1])*geo.nbVertices;
    int nbDofsPress = geo.nbVertices;
    int nbDofsVel = par.nbDofsPerNode()[0]*geo.nbVertices;
    int nbVertices = geo.nbVertices;

    /* local values */
    int nbDofsPerNode = par.nbDofsPerNode()[0]; 
    int nbNodesPerElement = geo.dimension()+1; 

    int m,n;
    code = MatGetOwnershipRange(C, &m, &n); CHKERR(code);

    Vec phi0_seq = getSequential(phi);

    if (world_rank == 0) cout << "NS: Assembling convection matrix." << endl;
    for (int partId = 0; partId < geo.elements().size(); partId++){
      for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on elements */

        vector<int> simplex = geo.elements()[partId][feId];
        /* get finite element coordinates */
        vector<vector<double>> coordinates(nbNodesPerElement);
        for (int i = 0; i < nbNodesPerElement; i++){
          coordinates[i] = geo.coordinates()[simplex[i]];
        }
        fe.setCoordinates(coordinates);
        fe.computeSize();

        if (feId % (geo.elements()[partId].size() / 5) == 0){
          if (world_rank == 0) cout << "  Elementary matrix for simplex: " << feId << "/" << geo.elements()[partId].size() - 1 << endl;
        }

        /* loc2glob dof mapping */
        vector<int> loc2glob(nbNodesPerElement*nbDofsPerNode);
        for (int i = 0; i < nbNodesPerElement; i++){
          for (int comp = 0; comp < nbDofsPerNode; comp++){
            loc2glob[nbDofsPerNode*i+comp] = nbDofsPerNode*simplex[i]+comp;
          }
        }

        /* Gather advection field on element */
        vector<double> u_el(nbNodesPerElement*nbDofsPerNode);
        code = VecGetValues(phi0_seq, nbNodesPerElement*nbDofsPerNode, &loc2glob[0], &u_el[0]); CHKERR(code);

        /* Set values from local to global */
        for (int j = 0; j < nbNodesPerElement; j++){
          vector<double> u_node_j(par.nbDofsPerNode()[0]);
          for (int i = 0; i < nbNodesPerElement; i++){
            for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
              u_node_j[comp] = u_el[par.nbDofsPerNode()[0]*j+comp];
            }
            for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
              if (par.nbDofsPerNode()[0]*simplex[i]+comp >= m && par.nbDofsPerNode()[0]*simplex[i]+comp < n){
                double Aij = par.density()*fe.advection(i, j, u_node_j);
                double SUPG = fe.size()*fe.size() * fe.supg(i, j, u_node_j);
                code = MatSetValue(C, par.nbDofsPerNode()[0]*simplex[i]+comp, par.nbDofsPerNode()[0]*simplex[j]+comp, Aij + SUPG, ADD_VALUES); CHKERR(code);
              }
            }
          }
        }
      }
    }
    code = MatAssemblyBegin(C, MAT_FINAL_ASSEMBLY); CHKERR(code);
    code = MatAssemblyEnd(C, MAT_FINAL_ASSEMBLY); CHKERR(code);
}
