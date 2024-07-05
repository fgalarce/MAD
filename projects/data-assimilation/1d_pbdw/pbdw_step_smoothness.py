from subprocess import *
from multiprocessing import Pool
from scipy import linalg
from tools import *
import matplotlib.pyplot as pl
import numpy as np
import fileinput
import scipy
import argparse
parser = argparse.ArgumentParser(description='parser')
parser.add_argument('-A', action="store", dest='A_gt', default=1)
parser.add_argument('-T', action="store", dest='T_gt', default=1)
args = parser.parse_args()
fontSize=16
fig_size_x=5
fig_size_y=4
#fig_size_x=10
#fig_size_y=8

# Parameters
nbSnapshots = 128
nbMeasures = 30
nbModes = 40
pixel_size = 0.04
im_start = 0.0 
im_end = 2.0*np.pi

# mesh
mesh = np.arange(0,2*np.pi, 0.01)
mesh_extended = np.arange(0,3*np.pi, 0.01)

# data assigned according to experiment
mani = 'step' # step, sinusoid
u_gt = np.zeros(len(mesh))
std_dev = 0
alpha = 0

def step(point, mesh, A = 1.0):
  heavy_side = np.zeros(len(mesh))
  for i in range(len(mesh)):
    if (mesh[i] > point):
      heavy_side[i] = A
  return heavy_side

if (mani == "step"):
  # manifold parameters
  Amplitude = np.random.uniform(0.8,1.2, nbSnapshots)
  step_size = np.random.uniform(1,2,nbSnapshots)
  T_start = 0.2*np.pi
  T_end = 0.4*np.pi
  delta_start = 0.0
  delta_end = 2.0*np.pi
  nbFrecuencies=10
  step_position = np.random.uniform(1, 3, nbSnapshots)

  # ground truth parameters
  A_gt = 1.0
  x_gt = 2.0
  T_gt = np.random.uniform(T_start, T_end, nbFrecuencies)
  std_dev = A_gt / 20 * 0
  alpha = A_gt/4 * 1
  delta_gt = np.random.uniform(delta_start, delta_end, nbFrecuencies)
  for j in range(nbFrecuencies):
    u_gt = u_gt + 1/5*np.sin(2*np.pi/T_gt[j]*mesh + delta_gt[j])/float(nbFrecuencies) 
  u_gt = u_gt + step(x_gt, mesh)
  manifold = np.zeros((len(mesh_extended), nbSnapshots))
  manifold_sin = np.zeros((len(mesh_extended), nbSnapshots))
  manifold_step = np.zeros((len(mesh_extended), nbSnapshots))
  for i in range(nbSnapshots):
    T = np.random.uniform(T_start, T_end, nbFrecuencies)
    delta = np.random.uniform(delta_start, delta_end, nbFrecuencies)
    for j in range(nbFrecuencies):
#      manifold_sin[:,i] = manifold_sin[:,i] + 1/8*np.sin(2*np.pi/T[j]*mesh_extended + delta[j])/float(nbFrecuencies)
      manifold_sin[:,i] = manifold_sin[:,i] + Amplitude[j]*np.sin(2*np.pi/T[j]*mesh_extended + delta[j])/float(nbFrecuencies)
    manifold_step[:,i] = step_size[i] * step(step_position[i], mesh_extended)
    manifold[:,i] = manifold_step[:,i] + manifold_sin[:,i]

# plot manifold
pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
for i in range(nbSnapshots):
  pl.plot(mesh, manifold_sin[0:len(mesh),i], alpha=0.7)
pl.xlabel("$x$", fontsize=fontSize)
pl.ylabel("snapshots", fontsize=fontSize)
show()

pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
for i in range(nbSnapshots):
  pl.plot(mesh, manifold_step[0:len(mesh),i], alpha=0.7)
pl.xlabel("$x$", fontsize=fontSize)
pl.ylabel("snapshots", fontsize=fontSize)
show()

pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
for i in range(nbSnapshots):
  pl.plot(mesh, manifold[0:len(mesh),i] + manifold_sin[0:len(mesh),i], alpha=0.7)
pl.xlabel("$x$", fontsize=fontSize)
pl.ylabel("snapshots", fontsize=fontSize)
show()


#
# Model reduction of smooth manifold
#
U_extended_smooth, Ssm, Vt = np.linalg.svd(manifold_sin, full_matrices=True, compute_uv=True)
U_extended_spike, Ssp, Vt = np.linalg.svd(manifold_step, full_matrices=True, compute_uv=True)
U_extended, S, Vt = np.linalg.svd(manifold, full_matrices=True, compute_uv=True)
U = np.zeros((len(mesh), nbModes))
Usm = np.zeros((len(mesh), nbModes))
Usp = np.zeros((len(mesh), nbModes))
pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
for i in range(nbModes):
  U[:,i] = U_extended[0:len(mesh),i]
  Usm[:,i] = U_extended_smooth[0:len(mesh),i]
  Usp[:,i] = U_extended_spike[0:len(mesh),i]
  pl.plot(mesh, U[:,i], label='mode ' + str(i+1))
pl.xlabel("$x$", fontsize=fontSize)
pl.ylabel("modes", fontsize=fontSize)
show(0)
print("is POD orthogonal: ", is_ortho(U))
print("is POD normalized: ", is_normal(U))

# plot model quality
pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
errorPOD = np.zeros(len(S))
errorPODsm = np.zeros(len(Ssm))
errorPODsp = np.zeros(len(Ssp))
for i in range(len(S)-1):
  errorPOD[i] = np.sqrt(sum(S[i+1:len(S)])) / np.sqrt(sum(S))
  errorPODsm[i] = np.sqrt(sum(Ssm[i+1:len(Ssm)])) / np.sqrt(sum(Ssm))
  errorPODsp[i] = np.sqrt(sum(Ssp[i+1:len(Ssp)])) / np.sqrt(sum(Ssp))
pl.plot(range(1,nbModes-2), errorPOD[0:nbModes-3])
pl.plot(range(1,nbModes-2), errorPODsm[0:nbModes-3])
pl.plot(range(1,nbModes-2), errorPODsp[0:nbModes-3])
pl.xlabel("$n$", fontsize=fontSize)
pl.ylabel("$\\epsilon_n$", fontsize=fontSize)
pl.yscale("log")
show()

# 
# Compute Riesz representers
#
dx = (im_end - im_start) / nbMeasures
rr = np.zeros((len(mesh), nbMeasures))
for i in range(nbMeasures):
  for j in range(len(mesh)):
#    if (mesh[j] >= im_start + i * pixel_size and mesh[j] < im_start + (i+1) * pixel_size):
    if (mesh[j] >= i * dx + im_start and mesh[j] < pixel_size + dx * i + im_start and mesh[j] < im_end):
      rr[j, i] = 1.0 
  rr[:,i] = rr[:,i] / np.linalg.norm(rr[:,i])
print("are rr orthogonal: ", is_ortho(rr))
print("are rr normalized: ", is_normal(rr))

# compute a-priori error bound
G = np.matmul(np.transpose(rr), U[:,0:nbSnapshots])
LHS = np.matmul(np.transpose(G), G)
UU, SS, VVt = np.linalg.svd(LHS, full_matrices=True, compute_uv=True)
pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
pl.plot(range(1,nbModes), SS[0:nbModes-1], label="beta_full")
# compute a-priori error bound
G = np.matmul(np.transpose(rr), Usm[:,0:nbSnapshots])
LHS = np.matmul(np.transpose(G), G)
UU, SS, VVt = np.linalg.svd(LHS, full_matrices=True, compute_uv=True)
pl.plot(range(1,nbModes), SS[0:nbModes-1], label="beta_sm")
pl.xlabel("$n$", fontsize=fontSize)
pl.ylabel("$\\beta$", fontsize=fontSize)
pl.yscale("log")
show(1)

#pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
#pl.plot(range(1,nbModes+1), errorPOD[0:nbModes]/SS[0:nbModes])
#pl.xlabel("$n$", fontsize=fontSize)
#pl.ylabel("$\\beta(V_n,W_m)^{-1} \\epsilon_n$", fontsize=fontSize)
#pl.yscale("log")
#show()
