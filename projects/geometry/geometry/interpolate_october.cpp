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

#define EPSILON_LOCAL 0.08 // fixed to current mesh size!

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
  cout << "Assembling stiffness for template geometry." << endl;
  Mat stiff_template = mat(nbDofs_i, nbDofs_i);
  MasterTetrahedron tet_template;
  tet_template.initialize(par);
  for (vector<int> tetra : geo_i.tetrahedron()[0]){ /* loop on tetra */
    /* get finite element coordinates */
    vector<vector<double>> coordinates(4);
    coordinates[0] = geo_i.coordinates()[tetra[0]];
    coordinates[1] = geo_i.coordinates()[tetra[1]];
    coordinates[2] = geo_i.coordinates()[tetra[2]];
    coordinates[3] = geo_i.coordinates()[tetra[3]];
    tet_template.setCoordinates(coordinates);
    /* Assemble elementary matrices */
    for (int j = 0; j < 4; j++){
      for (int i = 0; i < 4; i++){
        double Kij = tet_template.stiffness(i,j) ;
        for (int comp = 0; comp < 3; comp++){
          code = MatSetValue(stiff_template, 3*tetra[i]+comp, 3*tetra[j]+comp, Kij, ADD_VALUES); CHKERR(code);
        }
      }
    }
  }
  MatAssemblyBegin(stiff_template, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(stiff_template, MAT_FINAL_ASSEMBLY);

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
  cout << "Assembling stiffness for target geometry." << endl;
  Mat stiff = mat(nbDofs_j, nbDofs_j);
  MasterTetrahedron tet_target;
  tet_target.initialize(par);
  for (vector<int> tetra : geo_j.tetrahedron()[0]){ /* loop on tetra */
    /* get finite element coordinates */
    vector<vector<double>> coordinates(4);
    coordinates[0] = geo_j.coordinates()[tetra[0]];
    coordinates[1] = geo_j.coordinates()[tetra[1]];
    coordinates[2] = geo_j.coordinates()[tetra[2]];
    coordinates[3] = geo_j.coordinates()[tetra[3]];
    tet_target.setCoordinates(coordinates);
    /* Assemble elementary matrices */
    for (int j = 0; j < 4; j++){
      for (int i = 0; i < 4; i++){
        double Kij = tet_target.stiffness(i,j) ;
        for (int comp = 0; comp < 3; comp++){
          code = MatSetValue(stiff, 3*tetra[i]+comp, 3*tetra[j]+comp, Kij, ADD_VALUES); CHKERR(code);
        }
      }
    }
  }
  MatAssemblyBegin(stiff, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(stiff, MAT_FINAL_ASSEMBLY);

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

  MasterTriangle tri;
  tri.initialize(par);

  /*  - - - - - - - - - - - - - - - - - - - - -
        Compute I pod (\Omega_i)
      - - - - - - - - - - - - - - - - - - - - -  */
  cout << "\nComputing I pod(\Omega_i)" << endl;
  for (int k=0; k<Tgeo_j.coordinates().size(); k++){
    
    if (k % (Tgeo_j.coordinates().size() / 10) == 0){
      cout << "  Vertex " << k << " / " <<  Tgeo_j.coordinates().size() << endl;
    }

    bool found = false;
    for (vector<int> tetra : geo_i.tetrahedron()[0]){

      /* Set finite element */
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
              I_Tpod[idMode][3*k + comp] += pod[idMode][3*tetra[idNode] + comp] * bc[idNode];
            }
          }
        }
        found = true;
        break;
      }
    }
//    if (!found){
//      vector<double> P0 = Tgeo_j.coordinates(k);
//      
//      for (int partId = 0; partId < geo_i.triangles().size(); partId++) { /* loop on parts */
//        for (vector<int> tria : geo_i.triangles()[partId]){ /* loop on triangles */
//
//          /* get finite element coordinates */
//          vector<vector<double>> coordinates(3);
//          coordinates[0] = geo_i.coordinates()[tria[0]];
//          coordinates[1] = geo_i.coordinates()[tria[1]];
//          coordinates[2] = geo_i.coordinates()[tria[2]];
//
//          tri.setCoordinates(coordinates);
//          tri.computeNormal();
//
//          /* Compute projection of point in triangle */
//          vector<double> V1 = coordinates[0];
//          double cosalpha = dot(P0 - V1, tri.normal() ) / norm(P0 - V1);
//          vector<double> P0_projection = P0 - norm(V1 - P0) * cosalpha * tri.normal();
//
//          vector<int> neighbourd = geo_i.getTetraFromTriangle(tria);
//          vector<int> te = geo_i.tetrahedron()[neighbourd[0]][neighbourd[1]];
//
//          vector<vector<double>> coord(4);
//          coord[0] = geo_i.coordinates()[te[0]];
//          coord[1] = geo_i.coordinates()[te[1]];
//          coord[2] = geo_i.coordinates()[te[2]];
//          coord[3] = geo_i.coordinates()[te[3]];
//
//          tet.setCoordinates(coord);
//          vector<double> bc = tet.barycentric_coor(P0_projection);
//
//          if ( bc[0] >= -EPSILON_LOCAL && bc[1] >= -EPSILON_LOCAL && bc[2] >= -EPSILON_LOCAL && bc[3] >= -EPSILON_LOCAL) {
//            vector<double> bc_tetra = tet.barycentric_coor(P0);
//            for (int idMode = 0; idMode < par.nbModes(); idMode++){
//              for (int idNode = 0; idNode < 4; idNode++){
//                for (int comp = 0; comp < 3; comp++){
//                  I_Tpod[idMode][3*k + comp] += pod[idMode][3*te[idNode] + comp] * bc_tetra[idNode];
//                }
//              }
//            }
//            found = true;
//            break;
//          }
//        }
//      }
//    }

    if (!found){
      int nn = geo_i.nearest_neighbourd(Tgeo_j.coordinates(k));
      for (int idMode = 0; idMode < par.nbModes(); idMode++){
        for (int comp = 0; comp < 3; comp++){
          I_Tpod[idMode][3*k + comp] = pod[idMode][3*nn + comp];
        }
      }
      found = true;
    }
  }

  /* - - - - - - - - - - - - - - - 
      Compute I Tij
   - - - - - - - - - - - - - - - */
  cout << "\nComputing I Tij" << endl;
  for (int k=0; k<Tgeo_j.coordinates().size(); k++){
    
    if (k % (Tgeo_j.coordinates().size() / 10) == 0){
      cout << "  Vertex " << k << " / " <<  Tgeo_j.coordinates().size() << endl;
    }
    bool found = false;
    for (vector<int> tetra : geo_i.tetrahedron()[0]){

      /* Set finite element */
      vector<vector<double>> coord(4);
      coord[0] = geo_i.coordinates()[tetra[0]];
      coord[1] = geo_i.coordinates()[tetra[1]];
      coord[2] = geo_i.coordinates()[tetra[2]];
      coord[3] = geo_i.coordinates()[tetra[3]];

      tet.setCoordinates(coord);
      vector<double> bc = tet.barycentric_coor(Tgeo_j.coordinates()[k]);;

      if ( bc[0] >= 0 && bc[1] >= 0 && bc[2] >= 0 && bc[3] >= 0) {
        for (int idNode = 0; idNode < 4; idNode++){
          for (int comp = 0; comp < 3; comp++){
            I_Tij[3*k + comp] += Tij[3*tetra[idNode] + comp] * bc[idNode];
          }
        }
        found=true;
        break;
      }
    }
//    if (!found){
//      vector<double> P0 = Tgeo_j.coordinates(k);
//      
//      for (int partId = 0; partId < geo_i.triangles().size(); partId++) { /* loop on parts */
//        for (vector<int> tria : geo_i.triangles()[partId]){ /* loop on triangles */
//
//          /* get finite element coordinates */
//          vector<vector<double>> coordinates(3);
//          coordinates[0] = geo_i.coordinates()[tria[0]];
//          coordinates[1] = geo_i.coordinates()[tria[1]];
//          coordinates[2] = geo_i.coordinates()[tria[2]];
//
//          tri.setCoordinates(coordinates);
//          tri.computeNormal();
//
//          /* Compute projection of point in triangle */
//          vector<double> V1 = coordinates[0];
//          double cosalpha = dot(P0 - V1, tri.normal() ) / norm(P0 - V1);
//          vector<double> P0_projection = P0 - norm(V1 - P0) * cosalpha * tri.normal();
//
//          /* Check if projection lies in corresponding triangle */
//          vector<int> neighbourd = geo_i.getTetraFromTriangle(tria);
//          vector<int> te = geo_i.tetrahedron()[neighbourd[0]][neighbourd[1]];
//
//          vector<vector<double>> coord(4);
//          coord[0] = geo_i.coordinates()[te[0]];
//          coord[1] = geo_i.coordinates()[te[1]];
//          coord[2] = geo_i.coordinates()[te[2]];
//          coord[3] = geo_i.coordinates()[te[3]];
//
//          tet.setCoordinates(coord);
//          vector<double> bc = tet.barycentric_coor(P0_projection);
//
//          if ( bc[0] >= -EPSILON_LOCAL && bc[1] >= -EPSILON_LOCAL && bc[2] >= -EPSILON_LOCAL && bc[3] >= -EPSILON_LOCAL) {
//
//            vector<double> bc_tetra = tet.barycentric_coor(P0);
//
//            for (int idNode = 0; idNode < 4; idNode++){
//              for (int comp = 0; comp < 3; comp++){
//                I_Tij[3*k + comp] += Tij[3*te[idNode] + comp] * bc_tetra[idNode];
//              }
//            }
//            found = true;
//            break;
//          }
//        }
//      }
//    }
    if (!found){
      int nn = geo_i.nearest_neighbourd(Tgeo_j.coordinates(k));
      for (int comp = 0; comp < 3; comp++){
        I_Tij[3*k + comp] = Tij[3*nn + comp];
      }
      found = true;
    }
  }

  /* Write interpolated volumetric extension of LDDMM mapping */
  io_j.writeState(I_Tij, "I_Tij", 0.0);

  /* - - - - - - - - - - - - - - - - - - - - -
        Piola transform 
     - - - - - - - - - - - - - - - - - - - - - */
  vector<Vec> basis(par.nbModes());
  vector<Vec> basis_nopiola(par.nbModes());
  if (par.usePiola()){
    for (int idMode = 0; idMode < par.nbModes(); idMode++){
      vec(basis[idMode], nbDofs_j);
      basis[idMode] = cfd.piola(petsc(I_Tpod[idMode]), petsc(I_Tij));

      vec(basis_nopiola[idMode], nbDofs_j);
      code = VecCopy(petsc(I_Tpod[idMode]), basis_nopiola[idMode]); CHKERRQ(code);
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
    if (par.usePiola()){
      bc.block(basis_nopiola[idMode]);
    }
  }

  la.orthonormalize(basis);
  if (par.usePiola()){
    la.orthonormalize(basis_nopiola);
  }

  for (int idMode = 0; idMode < par.nbModes(); idMode++){
    if (par.usePiola()){
      io_j.writeState(basis_nopiola[idMode], "I_Tpod", idMode);
    }
    io_j.writeState(basis[idMode], "I_PTpod", idMode);
  }

  SlepcFinalize();
}
