#include <ultra-4d-flow.hpp>

/* function pointers are passed as argument to BC class */
vector<double> inlet(vector<double> x, double t){
  vector<double> bc(3, 0.0);
//  bc[2] = (1 - (x[0]*x[0] + x[1]*x[1])/(0.2*0.2))*sin(2*PI*t);
//  if (bc[2] < 0){
//    bc[2] = 0.0;  
//  }
  bc[2] = (1 - (x[0]*x[0] + x[1]*x[1])/(0.2*0.2));
  return bc;
}

//vector<double> inlet(vector<double> x, double t){
//  vector<double> bc(3, 0.0);
//  bc[2] = 1.0;
//  return bc;
//}

vector<double> noslip(vector<double> x, double t){
  vector<double> bc(3,0.0);
  return bc;
}

int main(int argc, char *argv[]){

  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;

  /* Parse data file */
  assert(argc == 2);
  string data_file = argv[1];
  Parameters par(data_file);
  par.print();

  /* Initialize MDMA objects */
  IO io;
  io.initialize(par);
    
  Geometry geo;
  geo.initialize(par, io);
  
  CFD cfd;
  cfd.initialize(par, geo, io);

  MasterTetrahedron tet; 
  tet.initialize(par);

  BoundaryConditions bc;
  bc.initialize(par, geo);

  int nbDofs = 4*io.nbVertices();
  int nbDofsPress = io.nbVertices();
  int nbDofsVel = 3*io.nbVertices();
  int nbVertices = io.nbVertices();

  Mat M = mat(nbDofs, nbDofs); /* Mass matrix */
  Mat B = mat(nbDofs, nbDofs); /* Static matrix */

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  int m,n;
  code = MatGetOwnershipRange(M, &m, &n); CHKERRQ(code);

  Tic tic;
  for (int tetraId = 0; tetraId < geo.tetrahedron()[0].size(); tetraId++){ /* loop on tetra */

    vector<int> tetra = geo.tetrahedron()[0][tetraId];
    /* get finite element coordinates */
    vector<vector<double>> coordinates(4);
    coordinates[0] = geo.coordinates()[tetra[0]];
    coordinates[1] = geo.coordinates()[tetra[1]];
    coordinates[2] = geo.coordinates()[tetra[2]];
    coordinates[3] = geo.coordinates()[tetra[3]];

    tet.setCoordinates(coordinates);
    tet.computeSize();

    if (tetraId % (geo.tetrahedron()[0].size() / 10) == 0){
      if (world_rank == 0) cout << "  \n* * * * * * * * * * * * * * * " << endl;
      if (world_rank == 0) cout << "  Elementary matrix for tetra: " << tetraId << "/" << geo.tetrahedron()[0].size() - 1 << endl;
      if (world_rank == 0) cout << "  Coordinates: " << endl << coordinates;
      if (world_rank == 0) cout << "  Labels: " << tetra << endl; 
      if (world_rank == 0) cout << "  Size: " << tet.size() << endl;
      if (world_rank == 0) cout << "  det(J) = " << tet.detJacobian() << endl; 
      if (world_rank == 0) cout << "  volume = " << tet.volume() << " cms^3." << endl;
      if (world_rank == 0) cout << "  Mass(1,1) = " << tet.mass(1,1) << endl;
      if (world_rank == 0) cout << "  Stiff(1,2) = " << tet.stiffness(1,2) << endl;
      if (world_rank == 0) cout << "  Mixed(2,2) = " << tet.mixed(2,2) << endl;
    }

    /* Assemble elementary matrices */
    for (int j = 0; j < 4; j++){
      for (int i = 0; i < 4; i++){
        for (int comp = 0; comp < 3; comp++){
          if (3*tetra[i]+comp >= m && 3*tetra[i]+comp < n){
            /* M */
            double Mij = par.density()/par.timeStep()*tet.mass(i,j);
            /* stiffness */
            double Kij = tet.stiffness(i,j) ;
            code = MatSetValue(M, 3*tetra[i]+comp, 3*tetra[j]+comp, Mij, ADD_VALUES); CHKERRQ(code);
            code = MatSetValue(B, 3*tetra[i]+comp, 3*tetra[j]+comp, par.viscosity() * Kij, ADD_VALUES); CHKERRQ(code);
            /* div_phi_i phi_j */
            double Bij = tet.mixed(i,j,comp);
            code = MatSetValue(B, 3*tetra[i]+comp, nbDofsVel + tetra[j], -1.0*Bij, ADD_VALUES); CHKERRQ(code); /* -B*/
          } 
          if (nbDofsVel + tetra[j] >= m && nbDofsVel + tetra[j] < n){
            /* div_phi_i phi_j */
            double Bij = tet.mixed(i,j,comp);
            code = MatSetValue(B, nbDofsVel + tetra[j], 3*tetra[i]+comp, 1.0*Bij, ADD_VALUES); CHKERRQ(code); /* B^T */
          }
        }
        if (nbDofsVel + tetra[i] >= m && nbDofsVel + tetra[i] < n){
          double Kij = tet.stiffness(i,j) ;
          /* Brezzi-Pitkaranta stabilization for P1-P1 elements */
          code = MatSetValue(B, nbDofsVel + tetra[i], nbDofsVel + tetra[j], tet.size()*tet.size() * Kij, ADD_VALUES); CHKERRQ(code);
        }
      }
    }
  }
  if (world_rank == 0) cout << endl << "Time elapsed: " << tic.toc() << " sec. "<< endl;
  code = MatAssemblyBegin(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(M, MAT_FINAL_ASSEMBLY); CHKERRQ(code);
  code = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY); CHKERRQ(code); 
  code = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY); CHKERRQ(code);

  if (world_rank == 0) cout << "Assembling mass matrix in the domain boundary." << endl;
  Mat Mb = mat(nbDofs, nbDofs);
  vector<double> normals_P0(3*geo.nbElementsBoundary(), 0.0);
  int idElement = 0;

  MasterTriangle tri;
  tri.initialize(par);

  /* Compute matrices only if they are not provided in the parameters file */
  cout << "CFD: Assembling mass matrix on the boundary. \n";
  for (int partId = 0; partId < geo.triangles().size(); partId++){ /* loop on parts */
    cout << "  Part: " << partId << endl;
    int triangleId = 0; 
    for (vector<int> tria : geo.triangles()[partId]){ /* loop on triangles */
      if ((triangleId % (geo.triangles()[partId].size() / 10)) == 0){
        cout << "  Triangle " << triangleId << " of " << geo.triangles()[partId].size() << endl;
      }
      /* get finite element coordinates */
      vector<vector<double>> coordinates(3);
      coordinates[0] = geo.coordinates()[tria[0]];
      coordinates[1] = geo.coordinates()[tria[1]];
      coordinates[2] = geo.coordinates()[tria[2]];

      tri.setCoordinates(coordinates);
      tri.computeJacobian();
      tri.computeElementMass();

      /* Here we store exclusively the normals for the outlets boundaries with a Windkessel model */
      if (geo.triangleLabels(partId) == par.outlet()){
        tri.computeNormal();
        for (int comp = 0; comp < 3; comp++){
          normals_P0[3*idElement + comp] = tri.normal()[comp];
        }
      }

      /* set values from local to global */
      for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
          double Mij = tri.mass(i,j);
          for (int comp = 0; comp < 3; comp++){
            code = MatSetValue(Mb, 3*tria[i]+comp, 3*tria[j]+comp, Mij, ADD_VALUES); CHKERRQ(code);
          }
        }
      }
      triangleId++;
      idElement++;
    }
  }
  cout << "CFD: Assembling global matrix" << endl;
  MatAssemblyBegin(Mb, MAT_FINAL_ASSEMBLY); 
  MatAssemblyEnd(Mb, MAT_FINAL_ASSEMBLY);

  Vec normals_P1 = vec(nbDofs);
  normals_P1 = cfd.interpolateP0_P1_boundary(normals_P0);

//  io.writeState(normals_P1, "normals");
//  exit(1);

  Vec phi0 = zeros(nbDofs);
  Vec phi = zeros(nbDofs);

  /* LHS matrix */
  Mat A = mat(nbDofs, nbDofs);
  code = MatDuplicate(M, MAT_COPY_VALUES, &A); CHKERRQ(code);
  code = MatAXPY(A, 1.0, B, DIFFERENT_NONZERO_PATTERN); CHKERRQ(code);

  KSP solver = configureKSP(A, "preonly", "lu", true);

  /* block LHS  */
  bc.time(0.0);
  bc.Dirichlet(par.inlet(), inlet);
  for (int i = 0; i < par.walls().size(); i++){
    bc.Dirichlet(par.walls()[i], noslip);
  }

  /* Save initial condition */
  io.writeState(phi, "phi", 0.0);

  /* block and assembly A */
  bc.block(A);

  ofstream file_out("./file.txt");

  double p_out;
  double p_d0 = 0.0;

  for (double t : range(0, par.nbIterations()*par.timeStep(), par.timeStep())){

    if (world_rank == 0) cout << "\n - - - - - - - - - - - - - - - - - - - - - - -  " << endl;
    if (world_rank == 0) cout << "  Solving stokes equation for time: " << t << endl;
    if (world_rank == 0) cout << "  Iteration: " << t/par.timeStep() << endl;

    /* Inertia from previous time step */
    Vec rhs = vec(nbDofs);
    MatMult(M, phi0, rhs);

    /* Apply tension boundary conditions with Windkessel model */
    double Q_out = cfd.flow(phi, par.outlet());
    cout << "  Coupling with 0D solver." << endl;
    double p_d = 1.0/(par.capacitance()/par.timeStep() + 1.0/par.distalResistance()) * (Q_out + p_d0*par.capacitance()/par.timeStep());
    p_out = p_d + par.resistance() * Q_out;
    p_d0 = p_d;
    cout << "  Flow = " << Q_out << "." << endl;
    cout << "  Distal pressure = " << p_d << "." << endl;
    cout << "  Proximal pressure = " << par.resistance() * Q_out << "." << endl;
    file_out << p_d << endl;

    Vec tension_out = zeros(nbDofs);
    for (int i = 0; i < nbDofsVel; i++){
      for (int comp = 0; comp < 3; comp++){
        double n_ij;
        int ind = 3*i+comp;
        VecGetValues(normals_P1, 1, &ind, &n_ij);
        VecSetValue(tension_out, ind, 1.0/p_out*n_ij, INSERT_VALUES);
      }
    }
    code = VecAssemblyBegin(tension_out); CHKERR(code);
    code = VecAssemblyEnd(tension_out); CHKERR(code);

    io.writeState(tension_out, "tension_out", t);
    Vec rhs_windkessel = vec(nbDofs);
    code = MatMult(Mb, tension_out, rhs_windkessel); CHKERR(code);
    code = VecAXPY(rhs, -1.0, rhs_windkessel); CHKERR(code);

    /* Apply Dirichlet boundary conditions */
    bc.time(t);
    bc.Dirichlet(par.inlet(), inlet);
    for (int i = 0; i < par.walls().size(); i++){
      bc.Dirichlet(par.walls()[i], noslip);
    }
    bc.block(rhs);

    if (world_rank == 0) cout << "  Solving linear system." << endl;
    code = KSPSolve(solver, rhs, phi); CHKERRQ(code); 

    double solNorm;
    VecNorm(phi, NORM_2, &solNorm);
    if (world_rank == 0) cout << "  l2-norm of solution: " << solNorm << endl;

    code = VecCopy(phi, phi0); CHKERRQ(code);


    io.writeState(phi, "phi", t);

    VecDestroy(&rhs);
//    VecDestroy(&tension_out);
  }

  par.finalize();
  file_out.close();
  KSPDestroy(&solver);
  MatDestroy(&A);
  VecDestroy(&phi);
  VecDestroy(&phi0);
  MatDestroy(&B);
  MatDestroy(&M);
  PetscFinalize();
}
