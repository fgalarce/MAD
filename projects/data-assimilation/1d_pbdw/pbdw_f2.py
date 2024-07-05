from subprocess import *
import matplotlib.pyplot as pl
import numpy as np
from multiprocessing import Pool
import fileinput
import scipy
from scipy import linalg

EPSILON = 1e-6

def getCenter(meshF2, fieldF2):
  x_start = 0.0
  i_start = 0
  for i in range(len(meshF2)):
    if (fieldF2[i] > 0):
      x_start = meshF2[i]
      i_start = i
      break

  x_end = 0.0
  for i in range(i_start, len(meshF2)):
    if (fieldF2[i] == 0):
      x_end = meshF2[i]
      i_end = i
      break;
  dirac = np.zeros(len(meshF2))
  dirac[i_start] = 1.0
  dirac[i_end] = 1.0
  center = x_start + (x_end - x_start)/2.0
#  pl.plot(meshF2, dirac)
  return center
    
def is_ortho(A):
  isortho = 1
  for i in range(len(A[0,:])):
    for j in range(len(A[0,:])):
      dp = np.dot(A[:,i], A[:,j])
      if (dp < - EPSILON and dp > EPSILON and i != j):
        isortho = 0
  return isortho

def is_normal(A):
  isnormal = 1
  for i in range(len(A[0,:])):
    dp = np.dot(A[:,i], A[:,i])
    if (dp < 1 - EPSILON and dp > 1 + EPSILON):
      isnormal = 0
  return isnormal

def wildcard(idIter):
  if (idIter < 10):
   iteration = "0000" + str(idIter);
  elif (idIter < 100):
    iteration = "000" + str(idIter);
  elif (idIter < 1000):
    iteration = "00" + str(idIter);
  elif (idIter < 10000):
    iteration = "0" + str(idIter);
  elif (idIter < 100000):
    iteration = str(idIter);
  return iteration

# Parameters
nbSnapshots = 64
nbModes = 1

r = 0.5 / np.cos(np.pi/4.0)
v0 = 1.0
alpha = 2
std_dev = v0 / 100 * 0.0

# mesh
dx = 0.002
mesh = np.arange(-2*r,2*r+dx, dx)

# synthetic measures parameters
pixel_size = 0.05
nbMeasures = 2*int(r/pixel_size)

# ground truth
u_gt = np.zeros(len(mesh)) 
for j in range(len(mesh)):
  if ((1-(mesh[j]/r)**2) > 0.0):
    u_gt[j] = v0*(1-(mesh[j]/r)**2)

# manifold generation
manifold = np.zeros((len(mesh), nbSnapshots))
pl.figure(num=None, figsize=(12, 8))
amplitude = np.arange(0.5,1.5,1.0/nbSnapshots)
for i in range(nbSnapshots):
  for j in range(len(mesh)):
    if (amplitude[i]*(1-(mesh[j]/r)**2) > 0):
      manifold[j,i] = amplitude[i]*(1-(mesh[j]/r)**2)
  pl.plot(mesh, manifold[:,i])
pl.tight_layout()
pl.grid()
pl.show()

# Model reduction
U, S, Vt = np.linalg.svd(manifold, full_matrices=True, compute_uv=True)
pl.figure(num=None, figsize=(12, 8))
for i in range(nbModes):
  pl.plot(U[:,i], label='mode ' + str(i+1))
pl.tight_layout()
pl.grid()
pl.title("POD modes")
pl.show()
#print("is POD orthogonal: ", is_ortho(U))
#print("is POD normalized: ", is_normal(U))

# plot model quality
pl.figure(num=None, figsize=(12, 8))

errorPOD = np.zeros(len(S))
for i in range(len(S)-1):
  errorPOD[i] = np.sqrt(sum(S[i+1:len(S)])) / np.sqrt(sum(S))
  
pl.plot(range(1,len(S)+1), errorPOD)
pl.yscale("log")
pl.title("Model quality")
pl.tight_layout()
pl.grid()
pl.show()

# Riesz representers
rr = np.zeros((len(mesh), nbMeasures))
pl.figure(num=None, figsize=(12, 8))
im_start = -r
im_end = im_start + pixel_size*nbMeasures
for i in range(nbMeasures):
  pointInside = 0
  for j in range(len(mesh)):
    if (mesh[j] >= im_start + i * pixel_size and mesh[j] < im_start + (i+1) * pixel_size):
      rr[j, i] = 1.0 
  rr[:,i] = rr[:,i] / np.linalg.norm(rr[:,i])
  pl.plot(mesh, rr[:,i])
pl.tight_layout()
pl.title("Riesz representers")
pl.grid()
pl.show()

print("are rr orthogonal: ", is_ortho(rr))
print("are rr normalized: ", is_normal(rr))

# synthetic measures
omega = np.matmul(np.transpose(rr), u_gt)
omega = np.matmul(rr, omega)
#for i in range(len(mesh)):
#  if (mesh[i] > im_start and mesh[i] < im_start + pixel_size*nbMeasures):
#    omega[i] = omega[i] + np.random.normal(-u_gt[i]/np.max(u_gt)*alpha, std_dev, 1)
pl.figure(num=None, figsize=(12, 8))
pl.plot(mesh, omega, label="$\omega$")
pl.tight_layout()
pl.grid()
pl.legend(loc='best')
pl.show()

# Load F2 Data
pl.figure(num=None, figsize=(12, 8))
theta = np.pi/4.0
F2 = np.loadtxt("../../../../data/measures/u_field2.txt")
meshF2 = np.loadtxt("../../../../data/measures/measuresF2_position.txt") 
nbLines = len(F2[0,:])
for i in range(nbLines):
  omegaF2 = F2[:,i] / np.cos(theta)
  meshF2 = meshF2 - getCenter(meshF2, omegaF2)
  pl.plot(meshF2, omegaF2, alpha=0.4)

pl.plot(mesh, u_gt, label="gt", linewidth=2.0)
pl.xlim([-2,2])
pl.plot(mesh, omega, label="$\omega$")
pl.tight_layout()
pl.grid()
pl.legend(loc='best')
pl.show()

omegaF2 = F2[:,0] / np.cos(theta)
voxelCenters = meshF2 - getCenter(meshF2, omegaF2)
# discard zero measures
index = 0
for i in range(len(omegaF2)):
  if (omegaF2[i] == 0):
    index = np.append(index, i)
voxelCenters = np.delete(voxelCenters, index)
omegaF2 = np.delete(omegaF2, index)

pl.plot(voxelCenters, omegaF2)

# Compute RR from data
pixel_size = meshF2[1] - meshF2[0]
print(pixel_size)
m = len(voxelCenters)
rrF2 = np.zeros((len(mesh),m))
for i in range(m):
  for j in range(len(mesh)):
    if (mesh[j] >= voxelCenters[i] - pixel_size/2.0 and mesh[j] < voxelCenters[i] + pixel_size/2.0):
      rrF2[j,i] = 1.0;
  rrF2[:,i] = rrF2[:,i] / np.linalg.norm(rrF2[:,i])
  pl.plot(mesh, rrF2[:,i])
pl.show()

print("are rrF2 orthogonal: ", is_ortho(rrF2))
print("are rrF2 normalized: ", is_normal(rrF2))

omega_synthetic = omega
omega = np.zeros(len(mesh))
for i in range(m):
  omega = omegaF2[i]*rrF2[:,i] + omega
rr=rrF2
nbMeasures=m

## compute a-priori error bound
#G = np.matmul(np.transpose(rr), U[:,0:nbSnapshots])
#LHS = np.matmul(np.transpose(G), G)
#UU, SS, VVt = np.linalg.svd(LHS, full_matrices=True, compute_uv=True)
#pl.figure(num=None, figsize=(12, 8))
## inf-suf constant
#pl.plot(range(1,nbMeasures), SS[0:nbMeasures-1])
#pl.yscale("log")
#pl.tight_layout()
#pl.grid()
#pl.show()
#
#pl.figure(num=None, figsize=(12, 8))
#pl.plot(range(1,nbMeasures-3), errorPOD[0:nbMeasures-4]/SS[0:nbMeasures-4])
#pl.tight_layout()
#pl.yscale("log")
#pl.grid()
#pl.show()

# normal equations
pl.figure(num=None, figsize=(12, 8))
u_star = np.zeros((len(mesh), nbModes))
u_star_synthetic = np.zeros((len(mesh), nbModes))
error_l2 = np.zeros(nbModes)
for j in range(nbModes):
  G = np.matmul(np.transpose(rr), U[:,0:j+1])
  LHS = np.matmul(np.transpose(G), G)
  print((U[:,0:j+1]).shape)
  print((omega.shape))
  RHS = np.matmul(np.transpose(U[:,0:j+1]), omega)
  RHS_synthetic = np.matmul(np.transpose(U[:,0:j+1]), omega_synthetic)
  c = np.linalg.solve(LHS, RHS)
  c_synthetic = np.linalg.solve(LHS, RHS_synthetic)
  for i in range(j+1):
    u_star[:,j] = u_star[:,j] + U[:,i] * c[i]
    u_star_synthetic[:,j] = u_star[:,j] + U[:,i] * c_synthetic[i]
  error_l2[j] = np.linalg.norm(u_star[:,j] - u_gt) / np.linalg.norm(u_gt)
#  pl.plot(mesh, u_star[:,j], label='u* (' + str(j) + ' modes)')
#pl.plot(mesh, u_gt, '*', label='ground truth', alpha=0.1)
#pl.legend(loc='best')
#pl.tight_layout()
#pl.grid()
#pl.show()

#pl.figure(num=None, figsize=(12, 8))
#pl.plot(range(nbModes), error_l2)
#pl.legend(loc='best')
#pl.yscale("log")
#pl.tight_layout()
#pl.grid()
#pl.show()

pl.figure(num=None, figsize=(12, 8))
pl.plot(mesh, u_gt, '*', label='ground truth', alpha=0.1)
pl.plot(mesh, u_star[:,np.argmin(error_l2)], label='u* (' + str(np.argmin(error_l2) +1) + ' modes)')
pl.plot(mesh, u_star_synthetic[:,np.argmin(error_l2)], label='u* (' + str(np.argmin(error_l2) +1) + ' modes) synthetic')
pl.legend(loc='best')
pl.tight_layout()
pl.grid()
pl.show()
