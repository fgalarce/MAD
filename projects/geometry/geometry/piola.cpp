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

#include<ultra-4d-flow.hpp>

int main(int argc, char *argv[]){

  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  /* Parse data file */
  assert(argc == 2);
  Parameters par(argv[1]);
  par.print();

  /* Initialize MDMA objects */
  IO io_i;
  io_i.initialize(par.templateGeometry());

  IO io_j;
  io_j.initialize(par);

  /* \Omega_i */
  Geometry geo_i;
  geo_i.initialize(io_i);     

  /* \Omega_j */
  Geometry geo_j;                 
  geo_j.initialize(par, io_j);

  Geometry Tgeo_i;   /* T(\Omega_i) */
  Geometry Tgeo_j;   /* T^(-1)\Omega_j */ 

  CFD cfd;
  cfd.initialize(par, geo_j, io_j);

  InnerProduct ip;
  ip.initialize(par, geo_j);

  LinearAlgebra la;
  la.initialize(par, ip);

  int nbVertices_i = io_i.nbVertices();
  int nbVertices_j = io_j.nbVertices();
  int nbDofs_i = 3*nbVertices_i;
  int nbDofs_j = 3*nbVertices_j; 

  /*
      Volumetric mapping \Omega_i --> \Omega_j
  */
  cout << "Harmonic extension of LDDMM template->target." << endl;
  Vec lddmm_ij = vec(nbDofs_i);
  loadVec(lddmm_ij, par.surfaceMapping());
  vector<int> bdNodes;
  for (int i = 0; i < geo_i.boundaryNodes().size(); i++){
    for (int j = 0; j < geo_i.boundaryNodes()[i].size(); j++){
      bdNodes.push_back( 3*geo_i.boundaryNodes()[i][j]+0 );
      bdNodes.push_back( 3*geo_i.boundaryNodes()[i][j]+1 );
      bdNodes.push_back( 3*geo_i.boundaryNodes()[i][j]+2 );
    }
  }
  Mat stiff_template = mat(nbDofs_i, nbDofs_i);
  loadMat(stiff_template, par.stiffness_matrix_path_template());
  blockMatrix(stiff_template, bdNodes);
  KSP ksp = configureKSP(stiff_template, "preonly", "lu");
  Vec x = vec(nbDofs_i);
  KSPSolve(ksp, lddmm_ij, x);
  KSPDestroy(&ksp);
  vector<double> Tij = stl(x); 
  vector<double> I_Tij(nbDofs_j, 0.0);

  /*
      Volumetric mapping \Omega_j --> \Omega_i
  */
  cout << "Harmonic extension of LDDMM target->template." << endl;
  Vec lddmm_ji = vec(nbDofs_j);
  loadVec(lddmm_ji, par.surfaceMappingInverse());
  bdNodes.clear();
  for (int i = 0; i < geo_j.boundaryNodes().size(); i++){
    for (int j = 0; j < geo_j.boundaryNodes()[i].size(); j++){
      bdNodes.push_back( 3*geo_j.boundaryNodes()[i][j] );
      bdNodes.push_back( 3*geo_j.boundaryNodes()[i][j]+1 );
      bdNodes.push_back( 3*geo_j.boundaryNodes()[i][j]+2 );
    }
  }
  Mat stiff = mat(nbDofs_j, nbDofs_j);
  loadMat(stiff, par.stiffness_matrix_path());
  blockMatrix(stiff, bdNodes);
  KSP ksp_inv = configureKSP(stiff, "preonly", "lu");
  Vec x_inv = vec(nbDofs_j);
  KSPSolve(ksp_inv, lddmm_ji, x_inv);
  KSPDestroy(&ksp_inv);
  vector<double> Tji = stl(x_inv);

  cout << "Initializing T(Omega_i) and T^{-1}(Omega_j)." << endl;
  Tgeo_i.initialize(io_i, Tij);
  Tgeo_j.initialize(io_j, Tji);

  cout << "Computing interpolation I(T_pod)." << endl;
  vector<vector<double>> pod(par.nbModes());
  vector<vector<double>> I_Tpod(par.nbModes());
  for (int idMode = 0; idMode < par.nbModes(); idMode++){
    pod[idMode].resize(nbDofs_i);
    I_Tpod[idMode].resize(nbDofs_j);
    fill(I_Tpod[idMode].begin(), I_Tpod[idMode].end(), 0.0);
    loadVec(pod[idMode], par.dirModel() + "/" + par.patientName() + "." + wildcard(idMode) + ".vct");
  }

  MasterTetrahedron tet;
  tet.initialize(par);

  for (int idMode = 0; idMode < par.nbModes(); idMode++){
    loadVec(I_Tpod[idMode], par.dirResults() + "/I_Tpod." + wildcard(idMode) + ".vct");
  }

  /* - - - - - - - - - - - - - - - 
      Compute I Tij
   - - - - - - - - - - - - - - - */
  for (int k=0; k<Tgeo_j.coordinates().size(); k++){
    
    if (k % (Tgeo_j.coordinates().size() / 10) == 0){
      cout << "  Vertex " << k << " / " <<  Tgeo_j.coordinates().size() << endl;
    }

    bool found = false; 
    for (int l=0; l<geo_i.tetrahedron()[0].size(); l++){

      /* Set finite element */
      vector<int> tetra = geo_i.tetrahedron()[0][l];
      vector<vector<double>> coord(4);
      coord[0] = geo_i.coordinates()[tetra[0]];
      coord[1] = geo_i.coordinates()[tetra[1]];
      coord[2] = geo_i.coordinates()[tetra[2]];
      coord[3] = geo_i.coordinates()[tetra[3]];

      tet.setCoordinates(coord);
      vector<double> bc = tet.barycentric_coor(Tgeo_j.coordinates()[k]);;

      if ( bc[0] >= 0 && bc[1] >= 0 && bc[2] >= 0 && bc[3] >= 0) {
        for (int idMode = 0; idMode < par.nbModes(); idMode++){
          for (int idNode = 0; idNode < 4; idNode++){
            for (int comp = 0; comp < 3; comp++){
              I_Tij[3*k + comp] += Tij[3*tetra[idNode] + comp] * bc[idNode];
            }
          }
        }
        found = true;
        break;
      }
    }
    if (!found){
      int nn = geo_i.nearest_neighbourd(Tgeo_j.coordinates(k));
      for (int idMode = 0; idMode < par.nbModes(); idMode++){
        for (int comp = 0; comp < 3; comp++){
          I_Tij[3*k + comp] = Tij[3*nn + comp];
        }
      }
      found = true;
    }
    if (!found){
      exit(1);
    }
  }
  io_j.writeState(I_Tij, "I_Tij");
  
//  loadVec(I_Tij, par.dirResults() + "/I_Tij.00000.vct");
  /* - - - - - - - - - - - - - - - - - - - - -
        Piola transform 
     - - - - - - - - - - - - - - - - - - - - - */
  vector<Vec> basis(par.nbModes());
  if (par.usePiola()){
    for (int idMode = 0; idMode < par.nbModes(); idMode++){
      vec(basis[idMode], nbDofs_j);
      basis[idMode] = cfd.piola(petsc(I_Tpod[idMode]), petsc(I_Tij));
    } 
  } else {
    for (int idMode = 0; idMode < par.nbModes(); idMode++){
      vec(basis[idMode], nbDofs_j);
      code = VecCopy(petsc(I_Tpod[idMode]), basis[idMode]); CHKERRQ(code);
    }
  }

  /* - - - - - - - - - - - - - - - - - - - - -
        Re orthonormalize and re-impose no-slip conditions 
    - - - - - - - - - - - - - - - - - - - - - */
  BoundaryConditions bc;
  bc.initialize(par, geo_j);
  bc.velocity(par.wall(), 0.0);

  for (int idMode = 0; idMode < par.nbModes(); idMode++){
    bc.block(basis[idMode]);
  }

  la.orthonormalize(basis);

  for (int idMode = 0; idMode < par.nbModes(); idMode++){
    io_j.writeState(basis[idMode], "piola", idMode);
  }

  SlepcFinalize();
}
