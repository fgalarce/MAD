from tools import *
 

def computeRiesRepresenters(mesh, nbMeasures, voxelCenters, pixel_size):
  rr = np.zeros((len(mesh), nbMeasures)) # normalized
  rr_nn = np.zeros((len(mesh), nbMeasures))
  for i in range(nbMeasures):
    pointInside = 0
    for j in range(len(mesh)):
      if (mesh[j] >= voxelCenters[i] - pixel_size/2.0 and mesh[j] < voxelCenters[i] + pixel_size/2.0):
        rr_nn[j, i] = 1.0 
    rr[:,i] = rr_nn[:,i] / np.linalg.norm(rr_nn[:,i])

#  print("are rr orthogonal: ", is_ortho(rr))
#  print("are rr normalized: ", is_normal(rr))
  return rr, rr_nn;

def assimilate(mesh, rr, U, omega, nbModes):
  u_star = np.zeros((len(mesh)))
  G = np.matmul(np.transpose(rr), U[:,0:nbModes])
  LHS = np.matmul(np.transpose(G), G)
  RHS = np.matmul(np.transpose(U[:,0:nbModes]), omega)
  c = np.linalg.solve(LHS, RHS)
  for i in range(nbModes):
    u_star = u_star + U[:,i] * c[i]
  return u_star

def assimilate_fingerprint(mesh, rr, U, omega, nbModes, manifold):
  nbSnapshots = len(manifold[0,:])
  # Noise filtering
  c_all = np.zeros((nbModes,nbSnapshots))
  c_box = np.zeros((nbModes,2))
  for i in range(nbSnapshots):
    c_all[:,i] = np.matmul(np.transpose(U[:,0:nbModes]), manifold[:,i])
  for i in range(nbModes):
    c_box[i,0] = np.min(c_all[i,:])
    c_box[i,1] = np.max(c_all[i,:])
  
  u_star = np.zeros((len(mesh)))
  G = np.matmul(np.transpose(rr), U[:,0:nbModes])
  LHS = np.matmul(np.transpose(G), G)
  RHS = np.matmul(np.transpose(U[:,0:nbModes]), omega)
  c = np.linalg.solve(LHS, RHS)
  for i in range(nbModes):
    # restrict coordinates within box
    if (c[i] < c_box[i,0]):
      c[i] = c_box[i,0]
    elif (c[i] > c_box[i,1]):
      c[i] = c_box[i,1]
    u_star = u_star + U[:,i] * c[i]
  return u_star

def data_mean(meshF2, experimental, rr, u):
  data_size = len(experimental[0,:])
  nbMeasures = len(experimental[:,0])
  theta = np.pi/4.0

  ratio_bias = 1.0 / max(u)
#  print(ratio_bias)
  experimental = experimental / ratio_bias

  data_mean = np.zeros(54)
  omegaF2coordData = np.zeros((nbMeasures,data_size))
  voxelCentersData = np.zeros((nbMeasures,data_size))
  for j in range(data_size):
    omegaF2coordData[:,j] = experimental[:,j] / np.cos(theta) 
    voxelCentersData[:,j] = meshF2 - getCenter(meshF2, omegaF2coordData[:,j])
    # discard zero measures
    index = 0
    for i in range(len(omegaF2coordData[:,j])):
      if (omegaF2coordData[i,j] == 0):
        index = np.append(index, i)
    voxelCenters = np.delete(voxelCentersData[:,j], index-2)
    omegaF2coord = np.delete(omegaF2coordData[:,j], index-2)
    data_mean = data_mean + 1.0/data_size*omegaF2coord[0:54]
    print(len(omegaF2coord))

  return data_mean
