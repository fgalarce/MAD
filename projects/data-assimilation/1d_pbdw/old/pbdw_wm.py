from subprocess import *
import matplotlib.pyplot as pl
import numpy as np
from multiprocessing import Pool
import fileinput
import scipy
from scipy import linalg
from tools import *
fontSize=16

# Parameters
nbSnapshots = 32
extra_modes = 1
nbMeasures = 14
nbModes = nbMeasures
pixel_size = 0.03
im_start = np.pi/4
im_end = 2*np.pi - np.pi/4
A_gt = 32.5
T_gt = 1.5*np.pi
std_dev = A_gt / 10 * 0
alpha = 10

# mesh
mesh = np.arange(0,3*np.pi, 0.01)

# ground truth
u_gt = A_gt * np.sin(2*np.pi/T_gt*mesh)

# manifold generation
A = np.arange(nbSnapshots)
T_start = 1.0*np.pi
T_end = 2.0*np.pi
T = np.arange(T_start, T_end, (T_end - T_start)/nbSnapshots)
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
pl.xlabel("x", fontsize=fontSize)
pl.ylabel("modes", fontsize=fontSize)
pl.grid()
pl.show()
print("is POD orthogonal: ", is_ortho(U))
print("is POD normalized: ", is_normal(U))

# noise collector
nc = np.zeros((len(mesh),nbSnapshots))
noise_iterations = 10
dx = (im_end - im_start) / nbMeasures
for i in range(nbSnapshots):
  for j in range(len(mesh)):
#    if (mesh[j] >= im_start  and mesh[j] < im_start + nbMeasures * pixel_size):
    if (mesh[j] >= i * dx + im_start and mesh[j] < pixel_size + dx * i + im_start and mesh[j] < im_end):
      for k in range(noise_iterations):
        if (np.max(manifold[:,i]) != 0):
          nc[j,i] = nc[j,i] + np.random.normal(-manifold[j,i]/np.max(manifold[:,i])*alpha, std_dev, 1)
  nc[:,i] = nc[:,i] / noise_iterations
eta = np.zeros((len(mesh)))
for i in range(1,np.shape(nc)[1]):
  eta = eta + nc[:,i]
eta = eta / nbSnapshots
pl.figure(num=None, figsize=(12, 8))
pl.plot(mesh, eta, label='nc mode ' + str(i+1))
pl.xlabel("x", fontsize=fontSize)
pl.ylabel("eta", fontsize=fontSize)
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
pl.xlabel("n", fontsize=fontSize)
pl.ylabel("error POD", fontsize=fontSize)
pl.yscale("log")
pl.tight_layout()
pl.grid()
pl.show()

# Riesz representers
rr = np.zeros((len(mesh), nbMeasures))
for i in range(nbMeasures):
  for j in range(len(mesh)):
#    if (mesh[j] >= im_start + i * pixel_size and mesh[j] < im_start + (i+1) * pixel_size):
    if (mesh[j] >= i * dx + im_start and mesh[j] < pixel_size + dx * i + im_start and mesh[j] < im_end):
      rr[j, i] = 1.0 
  rr[:,i] = rr[:,i] / np.linalg.norm(rr[:,i])
print("are rr orthogonal: ", is_ortho(rr))
print("are rr normalized: ", is_normal(rr))

omega = np.matmul(np.transpose(rr), u_gt)
omega = np.matmul(rr, omega)
for j in range(nbMeasures):
  for i in range(len(mesh)):
  #  if (mesh[i] > im_start and mesh[i] < im_start + pixel_size*nbMeasures):
    if (mesh[i] >= j * dx + im_start and mesh[i] < pixel_size + dx * j + im_start and mesh[i] < im_end):
      omega[i] = omega[i] + np.random.normal(-u_gt[i]/np.max(u_gt)*alpha, std_dev, 1)
pl.figure(num=None, figsize=(12, 8))
pl.plot(mesh, u_gt, label='GT')
pl.plot(mesh, omega, label='omega')
pl.xlabel("x", fontsize=fontSize)
pl.ylabel("solution", fontsize=fontSize)
pl.tight_layout()
pl.grid()
pl.show()

# correct RR with noise model 
for i in range(nbMeasures):
  eta = eta - proj(rr[:,i], eta)
eta = eta / np.sqrt(np.dot(eta, eta))
rr = np.column_stack((rr, eta))
print("are rr orthogonal: ", is_ortho(rr))
print("are rr normalized: ", is_normal(rr))
pl.figure(num=None, figsize=(12, 8))
for i in range(np.shape(rr)[1]):
  pl.plot(mesh, rr[:,i])
pl.xlabel("x", fontsize=fontSize)
pl.ylabel("RR + eta", fontsize=fontSize)
pl.tight_layout()
pl.grid()
pl.show()

fastplot(mesh, eta, "x", "eta")
fastplot(mesh, proj(eta, u_gt))

# compute a-priori error bound
G = np.matmul(np.transpose(rr), U[:,0:nbSnapshots])
LHS = np.matmul(np.transpose(G), G)
UU, SS, VVt = np.linalg.svd(LHS, full_matrices=True, compute_uv=True)
pl.figure(num=None, figsize=(12, 8))
# inf-suf constant
pl.plot(range(1,nbMeasures), SS[0:nbMeasures-1])
pl.xlabel("n", fontsize=fontSize)
pl.ylabel("beta", fontsize=fontSize)
pl.yscale("log")
pl.tight_layout()
pl.grid()
pl.show()

pl.figure(num=None, figsize=(12, 8))
pl.plot(range(1,nbMeasures-3), errorPOD[0:nbMeasures-4]/SS[0:nbMeasures-4])
pl.tight_layout()
pl.xlabel("n", fontsize=fontSize)
pl.ylabel("bound", fontsize=fontSize)
pl.yscale("log")
pl.grid()
pl.show()

# normal equations
pl.figure(num=None, figsize=(12, 8))
u_star = np.zeros((len(mesh), nbModes))
error_l2 = np.zeros(nbModes)
lamd = 1e-5
for j in range(1,nbModes):
  G = np.matmul(np.transpose(rr), U[:,0:j])
  LHS = np.matmul(np.transpose(G), G)
  LHS = LHS + lamd * np.identity(j)
  RHS = np.matmul(np.transpose(U[:,0:j]), omega)
  c = np.linalg.solve(LHS, RHS)
  for i in range(j):
    u_star[:,j] = u_star[:,j] + U[:,i] * c[i]
    error_l2[j-1] = np.linalg.norm(u_star[:,j] - u_gt) / np.linalg.norm(u_gt)
  pl.plot(mesh, u_star[:,j], label='u* (' + str(j) + ' modes)')
pl.plot(mesh, u_gt, '*', label='ground truth', alpha=0.1)
pl.plot(mesh, omega, label='omega', alpha=0.5)
pl.legend(loc='best')
pl.tight_layout()
pl.grid()
pl.xlabel("x", fontsize=fontSize)
pl.ylabel("solution", fontsize=fontSize)
pl.show()

pl.figure(num=None, figsize=(12, 8))
pl.plot(range(1,nbModes+1), error_l2)
pl.xlabel("n", fontsize=fontSize)
pl.ylabel("l2 error", fontsize=fontSize)
pl.yscale("log")
pl.tight_layout()
pl.grid()
pl.show()

Pu_Wm_eta = np.matmul(np.transpose(rr), u_gt)
Pu_Wm_eta = np.matmul(rr, Pu_Wm_eta)
pl.figure(num=None, figsize=(12, 8))
pl.plot(mesh, Pu_Wm_eta, label='Pu_Wm_eta')
pl.plot(mesh, omega, label='omega', alpha=0.5)
pl.plot(mesh, u_gt, '*', label='ground truth', alpha=0.1)
pl.plot(mesh, u_star[:,np.argmin(error_l2)], label='u* (' + str(np.argmin(error_l2)+1) + ' modes)')
pl.legend(loc='best')
pl.xlim(0, 2*np.pi)
pl.tight_layout()
pl.grid()
pl.show()
