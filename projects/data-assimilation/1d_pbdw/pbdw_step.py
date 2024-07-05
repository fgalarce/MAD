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

# Parameters
nbSnapshots = 2048
nbMeasures = 30
nbModes = 30
nbModes_omp = 4
pixel_size = 0.04
im_start = 0.0 
im_end = 2.0*np.pi

# mesh
mesh = np.arange(0,2*np.pi, 0.01)
mesh_extended = np.arange(0,3*np.pi, 0.01)

# data assigned according to experiment
mani = 'step' # step, sinusoid
utilde_method = 'omp' # least-squares, omp, neighbourd
u_gt = np.zeros(len(mesh))
std_dev = 0
alpha = 0

def step(point, mesh, A = 1.0):
  heavy_side = np.zeros(len(mesh))
  for i in range(len(mesh)):
    if (mesh[i] > point):
      heavy_side[i] = A
  return heavy_side

# manifold generation
if (mani == "sinusoid"):
  A_gt = 2.0
  T_gt = 0.82*np.pi
  std_dev = A_gt / 10 * 0
  alpha = A_gt/3 * 0
  #A_gt = args.A_gt
  #T_gt = args.T_gt
  u_gt = A_gt * np.sin(2*np.pi/T_gt*mesh)
  A_start = 1.0
  A_end = 10
  A = np.random.uniform(A_start, A_end, nbSnapshots)
  T_start = 0.75*np.pi
  T_end = 1.0*np.pi
  T = np.random.uniform(T_start, T_end, nbSnapshots)
  manifold = np.zeros((len(mesh_extended), nbSnapshots))
  for i in range(nbSnapshots):
    manifold[:,i] = A[i]*np.sin(2*np.pi/T[i]*mesh_extended)

if (mani == "step"):
  A_gt = 1.0
  T_start = 0.2*np.pi
  T_end = 0.4*np.pi
  delta_start = 0.0
  delta_end = 2.0*np.pi
  nbFrecuencies=10
  x_gt = 2.0
  T_gt = np.random.uniform(T_start, T_end, nbFrecuencies)
  std_dev = A_gt / 20 * 0
  alpha = A_gt/4 * 0
  step_position = np.random.uniform(1, 3, nbSnapshots)
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
      manifold_sin[:,i] = manifold_sin[:,i] + 1/8*np.sin(2*np.pi/T[j]*mesh_extended + delta[j])/float(nbFrecuencies)
    manifold_step[:,i] = step(step_position[i], mesh_extended)
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

# Model reduction 
U_extended, S, Vt = np.linalg.svd(manifold_sin, full_matrices=True, compute_uv=True)
U = np.zeros((len(mesh), nbModes))
pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
for i in range(nbModes):
  U[:,i] = U_extended[0:len(mesh),i]
  pl.plot(mesh, U[:,i], label='mode ' + str(i+1))
pl.xlabel("$x$", fontsize=fontSize)
pl.ylabel("modes", fontsize=fontSize)
show(1)
print("is POD orthogonal: ", is_ortho(U))
print("is POD normalized: ", is_normal(U))

# plot model quality
if (nbModes > 2):
  pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
  errorPOD = np.zeros(len(S))
  for i in range(len(S)-1):
    errorPOD[i] = np.sqrt(sum(S[i+1:len(S)])) / np.sqrt(sum(S))
  pl.plot(range(1,nbModes-2), errorPOD[0:nbModes-3])
  pl.xlabel("$n$", fontsize=fontSize)
  pl.ylabel("$\\epsilon_n$", fontsize=fontSize)
  pl.yscale("log")
  show()
  
# Riesz representers
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

omega = np.matmul(rr, np.matmul(np.transpose(rr), u_gt))
for j in range(nbMeasures):
  mean_gt = 0.0
  voxel_id = np.array([-1])
  for i in range(len(mesh)):
  #  if (mesh[i] > im_start and mesh[i] < im_start + pixel_size*nbMeasures):
    if (mesh[i] >= j * dx + im_start and mesh[i] < pixel_size + dx * j + im_start and mesh[i] < im_end):
      voxel_id = np.append(voxel_id, i)
      mean_gt = mean_gt + u_gt[i]
  mean_gt = mean_gt / (len(voxel_id)-1)
  omega[voxel_id[1:len(voxel_id)]] = omega[voxel_id[1:len(voxel_id)]] + np.random.normal(-mean_gt/np.max(u_gt)*alpha, std_dev, 1) * np.ones(len(voxel_id)-1)
pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
pl.plot(mesh, u_gt, label='$u_{true}$')
pl.plot(mesh, omega, label='$\\omega$')
pl.xlabel("$x$", fontsize=fontSize)
show(1)

#
# compute u tilde
#
u_tilde = np.zeros(len(mesh))
if (utilde_method == "neighbourd"):
  res = np.zeros(nbSnapshots)
  for i in range(nbSnapshots):
    res[i] = np.linalg.norm(omega - np.matmul(rr, np.transpose(rr)*manifold_step[0:len(mesh),i]))/np.linalg.norm(np.matmul(rr, np.transpose(rr)*manifold_step[0:len(mesh),i]))
  u_tilde = manifold_step[0:len(mesh),np.argmin(res)]

if (utilde_method == "least-squares"):
  pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
  v_star = np.zeros((len(mesh), nbModes))
  error_l2 = np.zeros(nbModes-1)
  for j in range(1,nbModes):
    G = np.matmul(np.transpose(rr), U[:,0:j])
    LHS = np.matmul(np.transpose(G), G)
    RHS = np.matmul(np.transpose(U[:,0:j]), omega)
    c = np.linalg.solve(LHS, RHS)
    for i in range(j):
      v_star[:,j] = v_star[:,j] + U[:,i] * c[i]
      error_l2[j-1] = np.linalg.norm(v_star[:,j] - u_gt) / np.linalg.norm(u_gt)
    pl.plot(mesh, v_star[:,j], label='u* (' + str(j) + ' modes)')
  u_tilde = v_star[0:len(mesh),np.argmin(error_l2)+1]
  show(1)

if (utilde_method == "omp"):
  pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
  v_omp = np.zeros((len(mesh), nbModes_omp))
  U_omp = np.zeros((len(mesh), nbModes_omp))
  z = np.zeros((len(mesh), nbSnapshots))
  error_l2 = np.zeros(nbModes_omp-1)
  # first function 
  for i in range(nbSnapshots):
    U_omp[:,0] = U_omp[:,0] + manifold[0:len(mesh),i]
  U_omp[:,0] = U_omp[:,0] / np.linalg.norm(U_omp[:,0])

#  nsnap = 0.0
#  id_first = 0
#  for i in range(nbSnapshots):
#    if (np.linalg.norm(manifold_step[:,i]) > nsnap):
#      id_first = i
#      nsnap = np.linalg.norm(manifold_step[:,i])
#  U_omp[:,0] = manifold_step[0:len(mesh),id_first]

  # measure manifold
  for i in range(nbSnapshots):
    measure = np.zeros(len(mesh))
    for j in range(nbMeasures):
      measure = measure + rr[:,j]*np.matmul(np.transpose(manifold_step[0:len(mesh),i]),rr[:,j])
    z[:,i] = measure / np.linalg.norm(measure)

  # run OMP 
  Pi_Wm_Ui = np.zeros((len(mesh), nbModes_omp))
  id_or = -1*np.ones(nbModes_omp-1).astype(int)
  for i in range(1,nbModes_omp):
    # build Pi_{W_m} U_i
    for j in range(i):
      for k in range(nbMeasures):
        Pi_Wm_Ui[:,j] = Pi_Wm_Ui[:,j] + rr[:,k]*np.matmul(np.transpose(U_omp[:,j]), rr[:,k])
    theta = np.zeros(len(mesh))
    for k in range(i):
      theta =  theta + Pi_Wm_Ui[:,k]*np.matmul(np.transpose(Pi_Wm_Ui[:,k]), omega)
    print(i)
    greedy = 0.0
    for j in range(nbSnapshots):
      res = 0.0
      res = np.abs(np.dot(omega - theta, z[:,j]))
      if (res > greedy):
        if ((isinside(j, id_or) == 0)):
          greedy = res
          id_or[i-1] = j
        print(j)
        print(res)
    U_omp[:,i] = manifold_step[0:len(mesh), id_or[i-1]] 
#    U_omp[:,0:i] = orthonormalize(U_omp[:,0:i])

  for i in range(nbModes_omp):
    pl.plot(mesh, U_omp[:,i], label="mode " + str(i))
#  pl.plot(mesh, U_omp[:,1], label="mode " + str(1))
#  pl.plot(mesh, U_omp[:,0], '*', label="mode " + str(0))
  print("------")
  print(id_or[nbModes_omp-2])
#  pl.plot(mesh, manifold_step[0:len(mesh),id_or[nbModes_omp-2]], '*', label="mode " + str(0))
  show(1)
#  u_tilde = U_omp[:,1]
#  u_tilde = manifold_step[0:len(mesh),]

if (utilde_method == "omp"):
  G = np.matmul(np.transpose(rr), U_omp[:,1])
  LHS = np.matmul(np.transpose(G), G)
  RHS = np.matmul(np.transpose(U_omp[:,1]), omega)
  c = RHS/LHS
  u_tilde = U_omp[:,1]*c
else: 
  # try reconstructions 
  for j in range(1,nbModes):
    G = np.matmul(np.transpose(rr), U_omp[:,0:j])
    LHS = np.matmul(np.transpose(G), G)
    RHS = np.matmul(np.transpose(U_omp[:,0:j]), omega)
    c = np.linalg.solve(LHS, RHS)
    for i in range(j):
      v_omp[:,j] = v_omp[:,j] + U_omp[:,i] * c[i]
      error_l2[j-1] = np.linalg.norm(v_omp[:,j] - u_gt) / np.linalg.norm(u_gt)
    pl.plot(mesh, v_omp[:,j], label='u* (' + str(j) + ' modes)')
  u_tilde = v_omp[:,np.argmin(error_l2)+1]

#
# eta
#
omega_eta = np.zeros(len(mesh))
for j in range(nbMeasures):
  mean_utilde = 0.0
  voxel_id = np.array([-1])
  for i in range(len(mesh)):
    if (mesh[i] >= j * dx + im_start and mesh[i] < pixel_size + dx * j + im_start and mesh[i] < im_end):
      voxel_id = np.append(voxel_id, i)
      mean_utilde = mean_utilde + u_tilde[i]
  mean_utilde = mean_utilde / (len(voxel_id)-1)
# expected value std_dev = 0
  omega_eta[voxel_id[1:len(voxel_id)]] = omega[voxel_id[1:len(voxel_id)]] + np.random.normal(mean_utilde/np.max(u_tilde)*alpha, 0.0, 1) * np.ones(len(voxel_id)-1)
pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
pl.plot(mesh, omega_eta, label='$\\omega + \\eta(\\omega)$')
pl.plot(mesh, omega, label='$\\omega$')
pl.plot(mesh, u_tilde, label='$\\tilde{u}$')
pl.plot(mesh, u_gt, label='$u_{true}$', alpha = 0.4)
pl.xlabel("$x$", fontsize=fontSize)
show(1)

pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
for i in range(np.shape(rr)[1]):
  pl.plot(mesh, rr[:,i])
pl.xlabel("$x$", fontsize=fontSize)
pl.ylabel("Riesz representers", fontsize=fontSize)
show()

# compute a-priori error bound
G = np.matmul(np.transpose(rr), U[:,0:nbSnapshots])
LHS = np.matmul(np.transpose(G), G)
UU, SS, VVt = np.linalg.svd(LHS, full_matrices=True, compute_uv=True)
pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
# inf-suf constant
pl.plot(range(1,nbModes), SS[0:nbModes-1])
pl.xlabel("$n$", fontsize=fontSize)
pl.ylabel("$\\beta$", fontsize=fontSize)
pl.yscale("log")
show(1)

pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
pl.plot(range(1,nbModes+1), errorPOD[0:nbModes]/SS[0:nbModes])
pl.xlabel("$n$", fontsize=fontSize)
pl.ylabel("$\\beta(V_n,W_m)^{-1} \\epsilon_n$", fontsize=fontSize)
pl.yscale("log")
show()

# compute shift corrector
zeta = np.zeros(len(mesh))
for i in range(nbMeasures):
  zeta = zeta + rr[:,i] * np.dot(rr[:,i], u_tilde)
pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
pl.plot(mesh, zeta, label='$\\zeta$')
show(1)

## normal equations
pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
u_star = np.zeros((len(mesh), nbModes))
error_l2 = np.zeros(nbModes-1)
for j in range(1,nbModes):
  G = np.matmul(np.transpose(rr), U[:,0:j])
  LHS = np.matmul(np.transpose(G), G)
  RHS = np.matmul(np.transpose(U[:,0:j]), omega_eta - zeta)
  c = np.linalg.solve(LHS, RHS)
  for i in range(j):
    u_star[:,j] = u_star[:,j] + U[:,i] * c[i]
    error_l2[j-1] = np.linalg.norm(u_star[:,j] - u_gt) / np.linalg.norm(u_gt)
  pl.plot(mesh, u_star[:,j], label='u* (' + str(j) + ' modes)')
pl.plot(mesh, u_gt, '*', label='$u_{true}$', alpha=0.1)
pl.plot(mesh, omega, label='$\\omega$', alpha=0.5)
pl.xlabel("x", fontsize=fontSize)
#pl.ylabel("solution", fontsize=fontSize)
pl.xlim(0, 2*np.pi)
pl.ylim(min(u_gt), max(u_gt))
show(1)

pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
pl.plot(range(1,nbModes), error_l2)
pl.xlabel("$n$", fontsize=fontSize)
pl.ylabel("$l^2$ relative error", fontsize=fontSize)
pl.yscale("log")
show()

pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
pl.plot(mesh, omega, label='$\\omega$', alpha=0.5)
pl.plot(mesh, u_gt, label='$u_{true}$', alpha=0.5)
pl.plot(mesh, u_tilde, label='$\\tilde{u}$')
pl.plot(mesh, u_star[:,np.argmin(error_l2)+1] + u_tilde, label='$u^*$ (' + str(np.argmin(error_l2)+1) + ' modes)')
#pl.plot(mesh, u_star[:,nbModes-1] + u_tilde, label='$u^*$ (' + str(nbModes) + ' modes)')
#pl.xlim(0, T)
pl.ylim(min(u_gt), max(u_gt))
show(1)
