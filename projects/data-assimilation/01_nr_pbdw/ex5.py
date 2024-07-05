from tools import *
from assimilation import *

# Model reduction
nbSnapshots = 64
nbModes = 4

# Ground truth
r = 0.5
v0 = 50.0

# mesh
dx = 0.001
x0 = 0.5
mesh = np.arange(-2*r+x0,2*r+dx+2*x0, dx)
n_gt = 1.0

u_gt = np.zeros(len(mesh)) 
for j in range(len(mesh)):
  if ((1-((mesh[j]-x0)/r)**2) > 0.0):
    u_gt[j] = v0*(1-np.fabs((mesh[j]-x0)/r)**(1+1/n_gt))

# -------------------------------------
# Manifold generation
manifold = np.zeros((len(mesh), nbSnapshots))
pl.figure(num=None, figsize=(4, 3))
amplitude = np.arange(40,60,(60-40)/nbSnapshots)
n_mani = np.random.uniform(0.8,1.2,nbSnapshots)
print(n_mani)
for i in range(nbSnapshots):
  for j in range(len(mesh)):
    velocity = amplitude[i]*(1-(np.fabs(mesh[j]-x0)/r)**(1+1/n_mani[i]))
    if (velocity > 0):
      manifold[j,i] = velocity
  pl.plot(mesh, manifold[:,i])
pl.tight_layout()
pl.grid("on")
pl.title("Manifold functions")
#pl.show()

# -------------------------------------
# Model reduction
U, S, Vt = np.linalg.svd(manifold, full_matrices=True, compute_uv=True)
pl.figure(num=None, figsize=(4, 3))
for i in range(nbModes):
  pl.plot(U[:,i], label='mode ' + str(i+1))
pl.tight_layout()
pl.grid("on")
#pl.show()
#print("is POD orthogonal: ", is_ortho(U))
#print("is POD normalized: ", is_normal(U))

# -------------------------------------
# Plot model quality
pl.figure(num=None, figsize=(4, 3))
errorPOD = np.zeros(len(S))
for i in range(len(S)-1):
  errorPOD[i] = np.sqrt(sum(S[i+1:len(S)])) / np.sqrt(sum(S))
pl.plot(range(1,len(S)+1), errorPOD)
pl.xlabel("n")
pl.ylabel("POD error")
pl.yscale("log")
pl.xlim([1,50])
pl.grid("on")
pl.tight_layout()
#pl.show()

# -------------------------------------
# Load LUC Data
pl.figure(num=None, figsize=(4, 3))
LUC = np.zeros(100)
ids = range(50,160)
error = np.zeros(len(ids))
for idMeasure in ids:
  fileMeasure = "/Users/fgalarce/home/01_research/MAD/data/measures/Doppler_LUC/Doppler_Data_28_03_24/doppler_luc."+ str(idMeasure) + ".txt"
  print("Reconstructing: " + fileMeasure)
  voxelCenters = np.arange(0.0, 1.0, 1.0/100)
  LUC_it = np.loadtxt(fileMeasure)[200:300]
  LUC = LUC_it
  # discard zero measures
  omegaCoord = LUC
  index = 0
  for i in range(len(omegaCoord)):
    if (omegaCoord[i] == 0):
      index = np.append(index, i)
  voxelCenters = np.delete(voxelCenters, index)
  omegaCoord = np.delete(omegaCoord, index)
  LUC = np.delete(LUC, index)
  
  # -------------------------------------
  # Riesz representers
  pixel_size = voxelCenters[1] - voxelCenters[0]
  nbMeasures = len(voxelCenters)
  rr, rr_nn = computeRiesRepresenters(mesh, nbMeasures, voxelCenters, pixel_size)
  omega = np.zeros(len(mesh))
  for i in range(nbMeasures):
    omega = rr_nn[:,i] * omegaCoord[i]  + omega # change F2 coordinates to mesh framework
  
  # -------------------------------------
  # Synthetic measures
  omega_synthetic = np.matmul(np.transpose(rr), u_gt)
  omega_synthetic = np.matmul(rr, omega_synthetic)
  
  # -------------------------------------
  # Compute reconstruction
  u_star_synt = assimilate(mesh, rr, U, omega_synthetic, nbModes) # synthetic
  u_tilde = assimilate_fingerprint(mesh, rr, U, omega, nbModes, manifold) # Field2
  
  # -------------------------------------
  # Measure u_tilde
  omega_tilde = np.matmul(np.transpose(rr), u_tilde)
  omega_tilde = np.matmul(rr, omega_tilde)
  
  # -------------------------------------
  # Compute simulated signal expected value
  F2 = np.loadtxt("../../../../data/measures/Doppler_F2/Field2_data.txt") # Load F2 Data
  meshF2 = np.arange(0.03, 1.03, 1.0/len(F2[:,0])) 
  nbLines = len(F2[0,:])
  for i in range(nbLines):
    omegaF2coord = F2[:,i]
  
  pixel_size_F2 = meshF2[1] - meshF2[0]
  rr_F2, rr_nn_f2 = computeRiesRepresenters(mesh, len(meshF2), meshF2, pixel_size_F2)
  omegaF2_average = np.zeros(len(F2[:,0]))
  measuresF2_average = np.zeros(len(mesh))
  for i in range(len(F2[0,:])):
    omegaF2_average = omegaF2_average + F2[:,i]*1.0/float(len(F2[0,:]))
  for i in range(len(rr_nn_f2[0,:])):
    measuresF2_average = measuresF2_average + rr_nn_f2[:,i]*omegaF2_average[i]
  omega_unbiased = omega_tilde + (omega_tilde - measuresF2_average)
  
  # -------------------------------------
  # Compute corrected reconstruction 
  u_star = assimilate_fingerprint(mesh, rr, U, omega_unbiased, nbModes, manifold) # synthetic
  print(np.linalg.norm(u_star-u_gt)/np.linalg.norm(u_gt))
  error[idMeasure-ids[0]] = np.linalg.norm(u_tilde-u_gt)/np.linalg.norm(u_gt)
np.savetxt("error_time_pbdw.txt", error)
