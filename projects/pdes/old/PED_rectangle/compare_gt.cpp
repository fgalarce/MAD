/*=============================================================================
  This file is part of the code MAD (or MAD)
  Multy-physics and Data assimilation for Medical Applications.
  Copyright (C) 2021,
    
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

#include <mad.hpp>

int main(int argc, char *argv[]){

  /* Deploy petsc and slepc */
  SlepcInitialize(&argc, &argv, NULL, NULL); 
  PetscErrorCode code;
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  /* Parse data file */
  assert(argc == 2);
  string data_file = argv[1]; 
  Parameters par(data_file);
  par.print();

  /* Initialize MAD objects */
  IO io;
  io.initialize(par);
    
  Geometry geo;
  geo.initialize(par, io);

  Vec u = vec(io.nbVertices()*par.nbDofsPerNode()[0]);
  loadVec(u, "./interpolated_solutions/u_184kv.00000.vct");

  Vec u_7 = vec(io.nbVertices()*par.nbDofsPerNode()[0]);
  loadVec(u_7, "./interpolated_solutions/Iu_7kv.00000.vct");
 
  Vec u_18 = vec(io.nbVertices()*par.nbDofsPerNode()[0]);
  loadVec(u_18, "./interpolated_solutions/Iu_18kv.00000.vct");
  
  Vec u_28 = vec(io.nbVertices()*par.nbDofsPerNode()[0]);
  loadVec(u_28, "./interpolated_solutions/Iu_28kv.00000.vct");

  Vec u_46 = vec(io.nbVertices()*par.nbDofsPerNode()[0]);
  loadVec(u_46, "./interpolated_solutions/Iu_46kv.00000.vct");

  VecAXPY(u_7, -1.0, u);
  VecAXPY(u_18, -1.0, u);
  VecAXPY(u_28, -1.0, u);
  VecAXPY(u_46, -1.0, u);

  Vec z = vec(io.nbVertices()*par.nbDofsPerNode()[0]);
  loadVec(z, "./interpolated_solutions/z_184kv.00000.vct");

  Vec z_7 = vec(io.nbVertices()*par.nbDofsPerNode()[0]);
  loadVec(z_7, "./interpolated_solutions/Iz_7kv.00000.vct");
 
  Vec z_18 = vec(io.nbVertices()*par.nbDofsPerNode()[0]);
  loadVec(z_18, "./interpolated_solutions/Iz_18kv.00000.vct");
  
  Vec z_28 = vec(io.nbVertices()*par.nbDofsPerNode()[0]);
  loadVec(z_28, "./interpolated_solutions/Iz_28kv.00000.vct");

  Vec z_46 = vec(io.nbVertices()*par.nbDofsPerNode()[0]);
  loadVec(z_46, "./interpolated_solutions/Iz_46kv.00000.vct");

  VecAXPY(z_7, -1.0, z);
  VecAXPY(z_18, -1.0, z);
  VecAXPY(z_28, -1.0, z);
  VecAXPY(z_46, -1.0, z);

  ofstream comparison_gt("./comparison_gt_u.txt");
  comparison_gt << norm(u_7);
  comparison_gt << endl;
  comparison_gt << norm(u_18);
  comparison_gt << endl;
  comparison_gt << norm(u_28);
  comparison_gt << endl;
  comparison_gt << norm(u_46);
  comparison_gt.close();

  Vec p = vec(io.nbVertices());
  loadVec(p, "./interpolated_solutions/p_184kv.00000.scl");

  Vec p_7 = vec(io.nbVertices());
  loadVec(p_7, "./interpolated_solutions/Ip_7kv.00000.scl");
 
  Vec p_18 = vec(io.nbVertices());
  loadVec(p_18, "./interpolated_solutions/Ip_18kv.00000.scl");
  
  Vec p_28 = vec(io.nbVertices());
  loadVec(p_28, "./interpolated_solutions/Ip_28kv.00000.scl");

  Vec p_46 = vec(io.nbVertices());
  loadVec(p_46, "./interpolated_solutions/Ip_46kv.00000.scl");

  VecAXPY(p_7, -1.0, p);
  VecAXPY(p_18, -1.0, p);
  VecAXPY(p_28, -1.0, p);
  VecAXPY(p_46, -1.0, p);

  ofstream comparison_gt_u("./comparison_gt_u.txt");
  comparison_gt_u << norm(u_7);
  comparison_gt_u << endl;
  comparison_gt_u << norm(u_18);
  comparison_gt_u << endl;
  comparison_gt_u << norm(u_28);
  comparison_gt_u << endl;
  comparison_gt_u << norm(u_46);
  comparison_gt_u.close();

  ofstream comparison_gt_z("./comparison_gt_z.txt");
  comparison_gt_z << norm(z_7);
  comparison_gt_z << endl;
  comparison_gt_z << norm(z_18);
  comparison_gt_z << endl;
  comparison_gt_z << norm(z_28);
  comparison_gt_z << endl;
  comparison_gt_z << norm(z_46);
  comparison_gt_z.close();

  ofstream comparison_gt_p("./comparison_gt_p.txt");
  comparison_gt_p << norm(p_7);
  comparison_gt_p << endl;
  comparison_gt_p << norm(p_18);
  comparison_gt_p << endl;
  comparison_gt_p << norm(p_28);
  comparison_gt_p << endl;
  comparison_gt_p << norm(p_46);
  comparison_gt_p.close();

  SlepcFinalize();
}
