from subprocess import *
import matplotlib.pyplot as pl
import numpy as np
from multiprocessing import Pool
import fileinput
import scipy
from scipy import linalg
from tools import *

# Parameters
nbSnapshots = 64
nbModes = 3
extra_modes = 1
nbMeasures = 30
pixel_size = 0.1
im_start = np.pi/4
im_end = im_start + pixel_size*nbMeasures
A_gt = 32.5
T_gt = 1.5*np.pi
std_dev = A_gt / 100
alpha = 2

# mesh
mesh = np.arange(0,3*np.pi, 0.01)

# ground truth
u_gt = A_gt * np.sin(2*np.pi/T_gt*mesh)

# manifold generation
A = np.arange(nbSnapshots)
T_start = 1.0*np.pi
T_end = 2.0*np.pi
T = np.arange(T_start, T_end, (T_end - T_start)/nbSnapshots)
print(np.shape(T))
manifold = np.zeros((len(mesh), nbSnapshots))
pl.figure(num=None, figsize=(12, 8))
for i in range(nbSnapshots):
  manifold[:,i] = A[i]*np.sin(2*np.pi/T[i]*mesh)
  pl.plot(mesh, manifold[:,i])
pl.tight_layout()
pl.grid()
pl.show()

# Model reduction
U, S, Vt = np.linalg.svd(manifold, full_matrices=True, compute_uv=True)
pl.figure(num=None, figsize=(12, 8))
for i in range(nbModes):
  pl.plot(U[:,i], label='mode ' + str(i+1))
pl.legend(loc='best')
pl.tight_layout()
pl.grid()
pl.show()
print("is POD orthogonal: ", is_ortho(U))
print("is POD normalized: ", is_normal(U))

# noise collector
nc = np.zeros((len(mesh),nbSnapshots*100))
for i in range(nbSnapshots*100):
  for j in range(len(mesh)):
    if (mesh[j] >= im_start  and mesh[j] < im_start + nbMeasures * pixel_size):
      nc[j,i] = np.random.normal(-u_gt[j]/np.max(u_gt)*alpha, std_dev, 1)
eta, eta_s, eta_Vt = np.linalg.svd(nc, full_matrices=True, compute_uv=True)
pl.figure(num=None, figsize=(12, 8))
for i in range(nbModes):
  pl.plot(eta[:,i], label='nc mode ' + str(i+1))
pl.legend(loc='best')
pl.tight_layout()
pl.grid()
pl.show()

# plot model quality
pl.figure(num=None, figsize=(12, 8))
errorPOD = np.zeros(len(S))
for i in range(len(S)-1):
  errorPOD[i] = np.sqrt(sum(S[i+1:len(S)])) / np.sqrt(sum(S))
pl.plot(range(1,len(S)+1), errorPOD)
pl.yscale("log")
pl.tight_layout()
pl.grid()
pl.show()

# Riesz representers
rr = np.zeros((len(mesh), nbMeasures))
pl.figure(num=None, figsize=(12, 8))
for i in range(nbMeasures):
  pointInside = 0
  for j in range(len(mesh)):
    if (mesh[j] >= im_start + i * pixel_size and mesh[j] < im_start + (i+1) * pixel_size):
      rr[j, i] = 1.0 
  rr[:,i] = rr[:,i] / np.linalg.norm(rr[:,i])
  pl.plot(mesh, rr[:,i])
pl.tight_layout()
pl.grid()
pl.show()
print("are rr orthogonal: ", is_ortho(rr))
print("are rr normalized: ", is_normal(rr))

# synthetic measures
omega = np.matmul(np.transpose(rr), u_gt)
omega = np.matmul(rr, omega)
for i in range(len(mesh)):
  if (mesh[i] > im_start and mesh[i] < im_start + pixel_size*nbMeasures):
    omega[i] = omega[i] + np.random.normal(-u_gt[i]/np.max(u_gt)*alpha, std_dev, 1)
pl.figure(num=None, figsize=(12, 8))
pl.plot(mesh, u_gt)
pl.plot(mesh, omega)
pl.tight_layout()
pl.grid()
pl.show()

# compute a-priori error bound
G = np.matmul(np.transpose(rr), U[:,0:nbSnapshots])
LHS = np.matmul(np.transpose(G), G)
UU, SS, VVt = np.linalg.svd(LHS, full_matrices=True, compute_uv=True)
pl.figure(num=None, figsize=(12, 8))
# inf-suf constant
pl.plot(range(1,nbMeasures), SS[0:nbMeasures-1])
pl.yscale("log")
pl.tight_layout()
pl.grid()
pl.show()

pl.figure(num=None, figsize=(12, 8))
pl.plot(range(1,nbMeasures-3), errorPOD[0:nbMeasures-4]/SS[0:nbMeasures-4])
pl.tight_layout()
pl.yscale("log")
pl.grid()
pl.show()

# normal equations
u_star = np.zeros((len(mesh)))
U = np.column_stack((U[:,0:nbModes],eta[:,0:extra_modes]))

for i in range(np.shape(U)[1]):
  pl.plot(U[:,i], label='mode ' + str(i+1))
pl.legend(loc='best')
pl.tight_layout()
pl.grid()
pl.show()

G = np.matmul(np.transpose(rr), U)
LHS = np.matmul(np.transpose(G), G)
RHS = np.matmul(np.transpose(U), omega)
c = np.linalg.solve(LHS, RHS)
for i in range(nbModes+extra_modes):
  u_star = u_star + U[:,i] * c[i]
#  error_l2[j-1] = np.linalg.norm(u_star[:,j] - u_gt) / np.linalg.norm(u_gt)
pl.plot(mesh, u_star, label='u* (' + str(j) + ' modes)')
pl.plot(mesh, u_gt, '*', label='ground truth', alpha=0.1)
pl.legend(loc='best')
pl.tight_layout()
pl.grid()
pl.show()

pl.figure(num=None, figsize=(12, 8))
pl.plot(mesh, u_gt, '*', label='ground truth', alpha=0.1)
pl.plot(mesh, u_star[:,np.argmin(error_l2)], label='u* (' + str(np.argmin(error_l2)) + ' modes)')
pl.legend(loc='best')
pl.tight_layout()
pl.grid()
pl.show()
