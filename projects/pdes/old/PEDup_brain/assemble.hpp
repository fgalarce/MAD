#include <mad.hpp>

void assembleSystem(Mat A, Mat  Mp, Mat  M_sol_1, Mat  M_sol_2, Parameters par, Geometry geo, MasterElement fe){

  PetscErrorCode code;
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  int nbDofsPerNode = par.nbDofsPerNode()[0]; 
  int nbVertices = geo.nbVertices;
  int nbDofsSol = nbVertices * nbDofsPerNode;
  int nbDofsPress = nbVertices;
  int nbDofs = nbDofsSol + nbDofsPress;

  Mat A = mat(nbDofs, nbDofs); /* Static matrix */

  /* local values */
  int nbNodesPerElement = geo.dimension()+1; 

  Vec mu_field = zeros(geo.nbTriangles());

  double Skemptom = 0.99;
  double BiotWillis = 1.0;

  int m,n;
  code = MatGetOwnershipRange(A, &m, &n); CHKERR(code);
  if (world_rank == 0) cout << "PORO: Assembling discretization matrix." << endl;
  for (int partId = 0; partId < geo.elements().size(); partId++){

    double lambda = par.youngModulus()[partId] * par.poissonRatio()[partId] / ((1.0 + par.poissonRatio()[partId])*(1.0 - 2.0*par.poissonRatio()[partId]));
    double mu = par.youngModulus()[partId]/(2.0*(1.0+par.poissonRatio()[partId]));
    double bulk  = par.youngModulus()[partId] / (3.0* ( 1.0 - 2.0*par.poissonRatio()[partId]));
    double c0 = (1.0 - BiotWillis * Skemptom) * BiotWillis / bulk / Skemptom; /* mass storage coefficient */
    double K_darcy =  par.permeability()[partId]/par.viscosity();

    if (world_rank == 0) cout << "  Assembling elementary matrices, part " << partId << endl;
    if (world_rank == 0) cout << "  mu = " << mu << endl;
    if (world_rank == 0) cout << "  lambda = " << lambda << endl;
    if (world_rank == 0) cout << "  c0 = " << c0 << endl;
    if (world_rank == 0) cout << "  Bulk = " << bulk << endl;
    if (world_rank == 0) cout << "  K / mu_f = " << K_darcy << endl;

    for (int feId = 0; feId < geo.elements()[partId].size(); feId++){ /* loop on elements */

      vector<int> simplex = geo.elements()[partId][feId];
      /* get finite element coordinates */
      vector<vector<double>> coordinates(nbNodesPerElement);
      for (int i = 0; i < nbNodesPerElement; i++){
        coordinates[i] = geo.coordinates()[simplex[i]];
      }

      fe.setCoordinates(coordinates);

      if (feId % (geo.elements()[partId].size() / par.verbose()) == 0){
        if (world_rank == 0) cout << "    Elementary matrix for element: " << feId << "/" << geo.elements()[partId].size() - 1 << endl;
      }

      /* Assemble elementary matrices */
      for (int j = 0; j < nbNodesPerElement; j++){
        for (int i = 0; i < nbNodesPerElement; i++){
          /* - - - - - -
            Newton law for solid 
          - - - - - - - */

          /* \epsilon(u) : \epsilon(v) (diagonal) */
          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
            if (nbDofsPerNode*simplex[i]+comp >= m && nbDofsPerNode*simplex[i]+comp < n){
              double Mij = par.density() * fe.mass(i,j) / (par.timeStep()*par.timeStep());
              double Bij = BiotWillis * fe.dphi_i_phi_j(j,i,comp);
              code = MatSetValue(M_sol_1, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp, Mij, ADD_VALUES); CHKERR(code);
              code = MatSetValue(A, nbDofsPerNode*simplex[i]+comp, nbDofsSol + simplex[j], Bij, ADD_VALUES); CHKERR(code); CHKERR(code);
            }
          }

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

          /* - - - - - -
            Darcy law for fluid  
          - - - - - - - */
          if (nbDofsSol + simplex[i] >= m && nbDofsSol + simplex[i] < n){
            double Kij = fe.stiffness(i,j);
            double Mij = fe.mass(i,j);
            code = MatSetValue(A, nbDofsSol + simplex[i], nbDofsSol + simplex[j], K_darcy * par.timeStep() * Kij, ADD_VALUES);
            code = MatSetValue(Mp, nbDofsSol + simplex[i], nbDofsSol + simplex[j], c0 * Mij, ADD_VALUES);
          }
          for (int comp = 0; comp < par.nbDofsPerNode()[0]; comp++){
            if (nbDofsSol + simplex[i] >= m && nbDofsSol + simplex[i] < n){
              double Bij = BiotWillis * fe.dphi_i_phi_j(j,i,comp);
              code = MatSetValue(M_sol_2, nbDofsSol + simplex[i], nbDofsPerNode*simplex[j]+comp, Bij , ADD_VALUES); CHKERR(code); /* - B^T */
            }
          }
        }
      }
    }
  }
  code = MatAssemblyBegin(M_sol_1, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(M_sol_1, MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyBegin(M_sol_2, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(M_sol_2, MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyBegin(Mp, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(Mp, MAT_FINAL_ASSEMBLY); CHKERR(code);
  code = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); CHKERR(code); 
  code = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); CHKERR(code);

  /* LHS static matrix */
  Mat C = mat(nbDofs, nbDofs);
  code = MatDuplicate(A, MAT_COPY_VALUES, &C); CHKERRQ(code);
  code = MatAXPY(C, 1.0, M_sol_1, DIFFERENT_NONZERO_PATTERN); CHKERRQ(code);
  code = MatAXPY(C, 1.0, M_sol_2, DIFFERENT_NONZERO_PATTERN); CHKERRQ(code);
  code = MatAXPY(C, 1.0, Mp, DIFFERENT_NONZERO_PATTERN); CHKERRQ(code);

}
