#include<mad.hpp>

Mat assembleMass(Parameters par, const Geometry & geo){

  PetscErrorCode code;

  int nbDofsPerNode = par.nbDofsPerNode()[0];
  int nbDofs = nbDofsPerNode*geo.nbVertices;
  int nbNodesPerElement = geo.dimension() + 1;
  Mat M = mat(nbDofs, nbDofs);

  MasterElement fe; 
  fe.initialize(par, geo.dimension());

  int low, high;
  code = MatGetOwnershipRange(M, &low, &high); CHKERR(code);
  for (int partId = 0; partId < geo.elements().size(); partId++){
    for (vector<int> simplex : geo.elements()[partId]){ /* loop on simplex */
      /* get finite element coordinates */
      vector<vector<double>> coordinates(geo.dimension()+1);
      for (int comp = 0; comp < geo.dimension() + 1; comp++){
        coordinates[comp] = geo.coordinates()[simplex[comp]];
      }
      fe.setCoordinates(coordinates);
      /* Assemble elementary matrices */
      for (int i = 0; i < nbNodesPerElement; i++){
        for (int j = 0; j < nbNodesPerElement; j++){
          for (int comp = 0; comp < nbDofsPerNode; comp++){
            if(nbDofsPerNode*simplex[i]+comp >= low && nbDofsPerNode*simplex[i]+comp < high){
              double Kij = fe.mass(i,j);
              code = MatSetValue(M, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp, Kij, ADD_VALUES); CHKERR(code);
            }
          }
        }
      }
    }
  }
  MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY);
  return M;
}

Mat assembleStiffness(Parameters par, const Geometry & geo){

  PetscErrorCode code;

  int nbDofsPerNode = par.nbDofsPerNode()[0];
  int nbDofs = nbDofsPerNode*geo.nbVertices;
  int nbNodesPerElement = geo.dimension() + 1;
  Mat K = mat(nbDofs, nbDofs);

  MasterElement fe; 
  fe.initialize(par, geo.dimension());

  int low, high;
  code = MatGetOwnershipRange(K, &low, &high); CHKERR(code);
  for (int partId = 0; partId < geo.elements().size(); partId++){
    for (vector<int> simplex : geo.elements()[partId]){ /* loop on simplex */
      /* get finite element coordinates */
      vector<vector<double>> coordinates(geo.dimension()+1);
      for (int comp = 0; comp < geo.dimension() + 1; comp++){
        coordinates[comp] = geo.coordinates()[simplex[comp]];
      }
      fe.setCoordinates(coordinates);
      /* Assemble elementary matrices */
      for (int i = 0; i < nbNodesPerElement; i++){
        for (int j = 0; j < nbNodesPerElement; j++){
          for (int comp = 0; comp < nbDofsPerNode; comp++){
            if(nbDofsPerNode*simplex[i]+comp >= low && nbDofsPerNode*simplex[i]+comp < high){
              double Kij = fe.stiffness(i,j);
              code = MatSetValue(K, nbDofsPerNode*simplex[i]+comp, nbDofsPerNode*simplex[j]+comp, Kij, ADD_VALUES); CHKERR(code);
            }
          }
        }
      }
    }
  }
  MatAssemblyBegin(K, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(K, MAT_FINAL_ASSEMBLY);
  return K;
}

vector<int> getBoundaryDofs(Parameters par, const Geometry & geo){
  vector<int> bdDofs;
  int nbDofsPerNode = par.nbDofsPerNode()[0];
  for (int i = 0; i < geo.boundaryNodes().size(); i++){
    for (int j = 0; j < geo.boundaryNodes()[i].size(); j++){
      for (int comp = 0; comp < nbDofsPerNode; comp++){
        bdDofs.push_back(nbDofsPerNode*geo.boundaryNodes()[i][j] + comp);
      }
    }
  }
  return bdDofs;
}

Vec harmonic_extension(Parameters par, const Geometry & geo, IO io, string disp, string mode = "direct"){

  int nbDofs = par.nbDofsPerNode()[0]*geo.nbVertices;
  Vec HE = zeros(nbDofs);

  Vec lddmm_ij = vec(nbDofs);
  io.loadVector(lddmm_ij, disp);
  vector<int> bdDofs = getBoundaryDofs(par, geo);
  Mat K = assembleStiffness(par, geo);

  /* scale to fit mesh */
  double scaling_factor;
  VecMax(lddmm_ij, NULL, &scaling_factor);
  VecScale(lddmm_ij, 1.0/scaling_factor);

  if (mode == "direct"){
    blockMatrix(K, bdDofs);
    KSP ksp = configureKSP(K, "gmres", "asm");
    KSPSolve(ksp, lddmm_ij, HE);
    KSPDestroy(&ksp);

  } else if (mode == "iterative"){
    Vec HE_0 = zeros(nbDofs);
    Mat M = assembleMass(par, geo);
    MatScale(K, par.viscosity());
    MatAXPY(K, 1.0, M, SAME_NONZERO_PATTERN);
    blockMatrix(K, bdDofs);
    KSP ksp;
    configureKSP(ksp, K, "gmres", "asm");
    for (int i = 0; i < par.nbIterations(); i++){
      MADprint("HE iteration ", i);
      Vec rhs = zeros(nbDofs);
      MatMult(M, HE_0, rhs);
      for (int j : bdDofs){
        vecSetInsert(rhs, j, stl(getSequential(lddmm_ij))[j]);
      }
      VecAssemblyBegin(rhs);
      VecAssemblyEnd(rhs);
      KSPSolve(ksp, rhs, HE);
      VecZeroEntries(HE_0);
      VecAXPY(HE_0, 1.0, HE);
    }
    KSPDestroy(&ksp);
  }

  VecScale(HE, scaling_factor);

  return HE;
}
