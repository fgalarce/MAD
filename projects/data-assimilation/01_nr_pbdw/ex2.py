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
fontSize=18
#fig_size_x=5
#fig_size_y=4
fig_size_x=5
fig_size_y=4

# Parameters
nbSnapshots = 1024
nbMeasures = 50
nbModes = 20
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

def compute_model_error(S):
  errorPOD = np.zeros(len(S))
  for i in range(len(S)-1):
    errorPOD[i] = np.sqrt(sum(S[i+1:len(S)])) / np.sqrt(sum(S))
  return errorPOD

def step(point, mesh, A = 1.0):
  heavy_side = np.zeros(len(mesh))
  for i in range(len(mesh)):
    if (mesh[i] > point):
      heavy_side[i] = A
  return heavy_side

if (mani == "step"):
  Amplitude = np.random.uniform(0.05,0.15, nbSnapshots)
  step_size = np.random.uniform(1,2,nbSnapshots)
  T_start = 0.2*np.pi
  T_end = 0.3*np.pi
  delta_start = 0.0
  delta_end = 2.0*np.pi
  nbFrecuencies=10
  A_gt = 0.1
  x_gt = 2.0
  T_gt = np.random.uniform(T_start, T_end, nbFrecuencies)
  std_dev = A_gt / 20 * 0
  alpha = 2*A_gt 
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
      manifold_sin[:,i] = manifold_sin[:,i] + Amplitude[j]*np.sin(2*np.pi/T[j]*mesh_extended + delta[j])/float(nbFrecuencies)
    manifold_step[:,i] = step_size[i] * step(step_position[i], mesh_extended)
    manifold[:,i] = manifold_step[:,i] + manifold_sin[:,i]

# plot manifold
pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
for i in range(10):
  pl.plot(mesh, manifold_sin[0:len(mesh),i], alpha=0.7)
pl.xlabel("$x$", fontsize=fontSize)
#pl.ylabel("snapshots", fontsize=fontSize)
pl.tight_layout()
pl.savefig("2023-07-30_ex2_manifold.pdf", format="pdf", bbox_inches="tight")
pl.show()

pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
for i in range(10):
  pl.plot(mesh, manifold_step[0:len(mesh),i], alpha=0.7)
pl.xlabel("$x$", fontsize=fontSize)
#pl.ylabel("snapshots", fontsize=fontSize)
pl.tight_layout()
pl.savefig("2023-07-30_ex2_manifold_slow.pdf", format="pdf", bbox_inches="tight")
pl.show()

pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
for i in range(10):
  pl.plot(mesh, manifold[0:len(mesh),i] + manifold_sin[0:len(mesh),i], alpha=0.7)
pl.xlabel("$x$", fontsize=fontSize)
#pl.ylabel("snapshots", fontsize=fontSize)
pl.tight_layout()
pl.savefig("2023-07-30_ex2_manifold_fast.pdf", format="pdf", bbox_inches="tight")
pl.show()

#
# Model reduction of full manifold
#
U_full_ext, S_full, Vt_full = np.linalg.svd(manifold, full_matrices=True, compute_uv=True)
U_full = np.zeros((len(mesh), nbModes))
pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
for i in range(nbModes):
  U_full[:,i] = U_full_ext[0:len(mesh),i]
for i in range(int(nbModes/2)):
  pl.plot(mesh, U_full[:,i], label='mode ' + str(i+1))
pl.xlabel("$x$", fontsize=fontSize)
#pl.ylabel("modes", fontsize=fontSize)
pl.savefig("2023-07-30_ex2_pod_fullManifold.pdf", format="pdf", bbox_inches="tight")
pl.tight_layout()
pl.show()

#
# Model reduction of fast manifold
#
U_extended, S, Vt = np.linalg.svd(manifold_sin, full_matrices=True, compute_uv=True)
U = np.zeros((len(mesh), nbModes))
pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
for i in range(nbModes):
  U[:,i] = U_extended[0:len(mesh),i]
for i in range(int(nbModes/2)):
  pl.plot(mesh, U[:,i], label='mode ' + str(i+1))
pl.xlabel("$x$", fontsize=fontSize)
#pl.ylabel("modes", fontsize=fontSize)
pl.savefig("2023-07-30_ex2_pod.pdf", format="pdf", bbox_inches="tight")
pl.tight_layout()
pl.show()
print("is POD orthogonal: ", is_ortho(U))
print("is POD normalized: ", is_normal(U))

#
# Model reduction of slow manifold
#
U_slow, S_slow, Vt_slow = np.linalg.svd(manifold_step, full_matrices=True, compute_uv=True)


pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
pl.plot(range(1,len(S)+1), compute_model_error(S), label="$\mathcal{M}_{fast}$")
pl.plot(range(1,len(S)+1), compute_model_error(S_full), label="$\mathcal{M}$", linewidth=1.5)
pl.plot(range(1,len(S)+1), compute_model_error(S_slow), "--", label="$\mathcal{M}_{slow}$")
pl.xlim(1, 40)
pl.ylim(1e-7, 1)
#pl.xlabel("$x$", fontsize=fontSize)
pl.tight_layout()
pl.yscale("log")
pl.legend(loc='best')
pl.tight_layout()
pl.show()

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
  pl.tight_layout()
  pl.savefig("2023-07-30_ex2_pod_fast_decay.pdf", format="pdf", bbox_inches="tight")
  pl.show()
 

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

pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
for i in range(np.shape(rr)[1]):
  pl.plot(mesh, rr[:,i])
pl.xlabel("$x$", fontsize=fontSize)
pl.ylabel("Riesz representers", fontsize=fontSize)
pl.tight_layout()
pl.savefig("2023-07-30_ex2_rr.pdf", format="pdf", bbox_inches="tight")
pl.show()

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
pl.legend(loc='best')
pl.savefig("2023-07-30_ex2_measures.pdf", format="pdf", bbox_inches="tight")
pl.tight_layout()
pl.show()


#
# compute flatener f
#
v_omp = np.zeros((len(mesh)))
U_omp = np.zeros(len(mesh))

# measure manifold
z = np.zeros((len(mesh), nbSnapshots))
for i in range(nbSnapshots):
  measure = np.zeros(len(mesh))
  for j in range(nbMeasures):
    measure = measure + rr[:,j]*np.matmul(np.transpose(manifold_step[0:len(mesh),i]),rr[:,j])
  z[:,i] = measure / np.linalg.norm(measure)

# run OMP 
Pi_Wm_Ui = np.zeros((len(mesh)))
# build Pi_{W_m} U_i
for k in range(nbMeasures):
  Pi_Wm_Ui = Pi_Wm_Ui + rr[:,k]*np.matmul(np.transpose(U_omp), rr[:,k])
greedy = 0.0
for j in range(nbSnapshots):
  res = 0.0
  res = np.abs(np.dot(omega, z[:,j]))
  if (res > greedy):
    greedy = res
    id_or = j
U_omp = manifold_step[0:len(mesh), id_or] 
pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
pl.xlabel("$x$", fontsize=fontSize)
pl.ylabel("$f_b$", fontsize=fontSize)
pl.plot(mesh, U_omp)
pl.tight_layout()
pl.show()

G = np.matmul(np.transpose(rr), U_omp)
LHS = np.matmul(np.transpose(G), G)
RHS = np.matmul(np.transpose(U_omp), omega)
c = RHS/LHS
f = U_omp*c

#
# compute u_tilde
#
pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
u_star = np.zeros((len(mesh), nbModes))
error_l2 = np.zeros(nbModes-1)
P_Wm_f = np.zeros(len(mesh))
for i in range(nbMeasures):
  P_Wm_f = P_Wm_f + np.dot(rr[:,i], f) * rr[:,i]
# normal equations
for j in range(1,nbModes):
  G = np.matmul(np.transpose(rr), U[:,0:j])
  LHS = np.matmul(np.transpose(G), G)
  RHS = np.matmul(np.transpose(U[:,0:j]), omega - P_Wm_f)
  c = np.linalg.solve(LHS, RHS)
  for i in range(j):
    u_star[:,j] = u_star[:,j] + U[:,i] * c[i]
    error_l2[j-1] = np.linalg.norm(u_star[:,j] - u_gt) / np.linalg.norm(u_gt)
  pl.plot(mesh, u_star[:,j], label='u* (' + str(j) + ' modes)')
#u_tilde = u_star[:,np.argmin(error_l2)+1] + f
u_tilde = u_star[:,nbModes-1] + f

pl.plot(mesh, u_gt, '*', label='$u_{true}$', alpha=0.2)
pl.plot(mesh, omega, label='$\\omega$', alpha=0.5)
pl.xlim(0, 2*np.pi)
pl.ylim(min(u_gt), max(u_gt))
pl.tight_layout()
pl.legend(loc='best')
pl.show()

pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
pl.plot(mesh, u_gt, '*', label='$u_{true}$', alpha=0.2)
pl.plot(mesh, omega, label='$\\omega$', alpha=0.5)
pl.plot(mesh, f, label='$f$', alpha=0.5)
pl.plot(mesh, u_star[:,nbModes-1] + f, label='$\\tilde{u}_f$ ' + str(np.argmin(error_l2)+1) + " modes" , alpha=0.5)
pl.plot(mesh, u_tilde, label='$\\tilde{u}_f$ ', alpha=0.5)
pl.xlabel("x", fontsize=fontSize)
pl.tight_layout()
pl.legend(loc='best')
pl.show()

pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
pl.plot(range(1,nbModes), error_l2)
pl.xlabel("$n$", fontsize=fontSize)
pl.ylabel("$l^2$ relative error", fontsize=fontSize)
pl.yscale("log")
pl.tight_layout()
pl.show()

# -------------------------------------
# eta
#
eta = np.zeros(len(mesh))
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
  eta[voxel_id[1:len(voxel_id)]] = eta[voxel_id[1:len(voxel_id)]] + np.random.normal(mean_utilde/np.max(f)*alpha, 0.0, 1) * np.ones(len(voxel_id)-1)
  omega_eta[voxel_id[1:len(voxel_id)]] = omega[voxel_id[1:len(voxel_id)]]
omega_eta = omega_eta + eta
pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
pl.plot(mesh, omega_eta, label='$\\eta(\\omega)$')
pl.plot(mesh, omega, label='$\\omega$')
pl.plot(mesh, f, label='$f^*$')
pl.plot(mesh, u_gt, label='$u_{gt}$', alpha = 0.4)
pl.xlabel("$x$", fontsize=fontSize)
pl.legend(loc='best')
pl.tight_layout()
pl.show()

# -------------------------------------
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
pl.tight_layout()
pl.show()

pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
pl.plot(range(1,nbModes+1), errorPOD[0:nbModes]/SS[0:nbModes])
pl.xlabel("$n$", fontsize=fontSize)
pl.ylabel("$\\beta(V_n,W_m)^{-1} \\epsilon_n$", fontsize=fontSize)
pl.yscale("log")
pl.tight_layout()
pl.show()


# -------------------------------------
#  Compute unbiased flattener f_u 
for j in range(nbSnapshots):
  res = 0.0
  res = np.abs(np.dot(omega + eta, z[:,j]))
  if (res > greedy):
    greedy = res
    id_or = j
U_omp = manifold_step[0:len(mesh), id_or] 
G = np.matmul(np.transpose(rr), U_omp)
LHS = np.matmul(np.transpose(G), G)
RHS = np.matmul(np.transpose(U_omp), omega + eta)
c = RHS/LHS
f_u = U_omp*c
pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
pl.plot(mesh, f, label="$f^*$", linewidth=3, alpha=0.8)
pl.plot(mesh, f_u, label="$f_u$")
pl.plot(mesh, u_gt, label='$u_{gt}$')
pl.plot(mesh, u_tilde - f, label='$\\tilde{u}_f$')
pl.plot(mesh, u_tilde - f + f_u, label='$u^*$ sPBDW')
pl.xlabel("$x$", fontsize=fontSize)
pl.legend(loc='best')
pl.tight_layout()
pl.show()

pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
u_star = np.zeros((len(mesh), nbModes))
error_l2 = np.zeros(nbModes-1)
zeta = np.zeros(len(mesh))
for i in range(nbMeasures):
  zeta = zeta + rr[:,i]*np.dot(rr[:,i], f_u)
# normal equations
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
#pl.xlabel("x", fontsize=fontSize)
#pl.ylabel("solution", fontsize=fontSize)
pl.xlim(0, 2*np.pi)
pl.ylim(min(u_gt), max(u_gt))
pl.tight_layout()
pl.show()

pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
pl.plot(range(1,nbModes), error_l2)
pl.xlabel("$n$", fontsize=fontSize)
pl.ylabel("$l^2$ relative error", fontsize=fontSize)
pl.yscale("log")
pl.tight_layout()
pl.show()

#
# compute u PBDW 
#
u_star_pbdw = np.zeros((len(mesh), nbModes))
error_l2_pbdw = np.zeros(nbModes-1)
P_Wm_f_pbdw = np.zeros(len(mesh))
# normal equations
for j in range(1,nbModes):
  G = np.matmul(np.transpose(rr), U_full[:,0:j])
  LHS = np.matmul(np.transpose(G), G)
  RHS = np.matmul(np.transpose(U[:,0:j]), omega)
  c = np.linalg.solve(LHS, RHS)
  for i in range(j):
    u_star_pbdw[:,j] = u_star_pbdw[:,j] + U_full[:,i] * c[i]
    error_l2_pbdw[j-1] = np.linalg.norm(u_star_pbdw[:,j] - u_gt) / np.linalg.norm(u_gt)


pl.figure(num=None, figsize=(fig_size_x, fig_size_y))
#pl.plot(mesh, omega, label='$\\omega$', alpha=0.5)
pl.plot(mesh, u_gt, "--", label='$u_{true}$')
#pl.plot(mesh, u_tilde, label='$\\tilde{u}$')
#pl.plot(mesh, u_star[:,np.argmin(error_l2)+1] + f_u, label='$u^*$ (' + str(np.argmin(error_l2)+1) + ' modes)')
pl.plot(mesh, u_star[:,nbModes-1] + f_u, label='sPBDW')
pl.plot(mesh, u_star_pbdw[:,nbModes-1], label='PBDW', alpha=0.2)
#pl.plot(mesh, u_tilde - f + f_u, label='$u^*$?')
#pl.plot(mesh, u_star[:,nbModes-1] + u_tilde, label='$u^*$ (' + str(nbModes) + ' modes)')
#pl.xlim(0, T)
pl.ylim([min(u_gt)-0.1, max(u_gt)+0.1])
pl.tight_layout()
pl.legend(loc='best')
pl.savefig("2023-07-30_ex2_comparison_wrt_pbdw.pdf", format="pdf", bbox_inches="tight")
pl.show()
