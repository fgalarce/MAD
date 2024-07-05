from tools import *
from assimilation import *

# -------------------------------------
# Parameters

# Model reduction
nbSnapshots = 64
nbModes = 1

# Ground truth
r = 0.5 / np.cos(np.pi/4.0)
v0 = 1.0

# mesh
dx = 0.002
mesh = np.arange(-2*r,2*r+dx, dx)

u_gt = np.zeros(len(mesh)) 
for j in range(len(mesh)):
  if ((1-(mesh[j]/r)**2) > 0.0):
    u_gt[j] = v0*(1-(mesh[j]/r)**2)

# -------------------------------------
# Manifold generation
manifold = np.zeros((len(mesh), nbSnapshots))
pl.figure(num=None, figsize=(6, 4))
amplitude = np.arange(0.5,1.5,1.0/nbSnapshots)
for i in range(nbSnapshots):
  for j in range(len(mesh)):
    if (amplitude[i]*(1-(mesh[j]/r)**2) > 0):
      manifold[j,i] = amplitude[i]*(1-(mesh[j]/r)**2)
  pl.plot(mesh, manifold[:,i])
pl.tight_layout()
pl.grid()
pl.show()

# -------------------------------------
# Model reduction
U, S, Vt = np.linalg.svd(manifold, full_matrices=True, compute_uv=True)
pl.figure(num=None, figsize=(6, 4))
for i in range(nbModes):
  pl.plot(U[:,i], label='mode ' + str(i+1))
pl.tight_layout()
pl.grid()
pl.title("POD modes")
pl.show()
#print("is POD orthogonal: ", is_ortho(U))
#print("is POD normalized: ", is_normal(U))

# -------------------------------------
# Plot model quality
pl.figure(num=None, figsize=(6, 4))
errorPOD = np.zeros(len(S))
for i in range(len(S)-1):
  errorPOD[i] = np.sqrt(sum(S[i+1:len(S)])) / np.sqrt(sum(S))
pl.plot(range(1,len(S)+1), errorPOD)
pl.yscale("log")
pl.title("Model quality")
pl.tight_layout()
pl.grid()
pl.show()

# -------------------------------------
# Load F2 Data
theta = np.pi/4.0
F2 = np.loadtxt("../../../../data/measures/Doppler_F2/u_field2.txt")
meshF2 = np.loadtxt("../../../../data/measures/Doppler_F2/measuresF2_position.txt") 
nbLines = len(F2[0,:])
for i in range(nbLines):
  omegaF2coord = F2[:,i] / np.cos(theta)
  meshF2 = meshF2 - getCenter(meshF2, omegaF2coord)
  pl.plot(meshF2, omegaF2coord, alpha=0.4)
pl.plot(mesh, u_gt, label="$u_{gt}$", linewidth=2.0)
pl.xlim([-1.5*r, 1.5*r])
pl.xlabel("x [cm]")
pl.ylabel("u [m/s]")
pl.grid("on")
pl.legend(loc='best')
pl.tight_layout()
pl.show()

omegaF2coord = F2[:,0] / np.cos(theta) 
voxelCenters = meshF2 - getCenter(meshF2, omegaF2coord)
# discard zero measures
index = 0
for i in range(len(omegaF2coord)):
  if (omegaF2coord[i] == 0):
    index = np.append(index, i)
voxelCenters = np.delete(voxelCenters, index)
omegaF2coord = np.delete(omegaF2coord, index)
limitMaxMeasures = 54 # limit to maximal number to makie it coherent with data-base as a whole
voxelCenters = voxelCenters[0:limitMaxMeasures]
omegaF2coord = omegaF2coord[0:limitMaxMeasures]

# -------------------------------------
# Riesz representers
pixel_size = meshF2[1] - meshF2[0]
nbMeasures = len(voxelCenters)
print("pixel size in F2 data: ", pixel_size)
rr, rr_nn = computeRiesRepresenters(mesh, nbMeasures, voxelCenters, pixel_size)
omegaF2 = np.zeros(len(mesh))
pl.figure(num=None, figsize=(6, 4))
for i in range(nbMeasures):
  omegaF2 = rr_nn[:,i] * omegaF2coord[i]  + omegaF2 # change F2 coordinates to mesh framework
  pl.plot(mesh, rr[:,i])
pl.plot(voxelCenters, omegaF2coord, "*", label="$\\omega$ coordinates")
pl.plot(mesh, omegaF2, label="$\\omega$")
pl.xlim([-1.5*r, 1.5*r])
pl.legend(loc='best')
pl.xlabel("x [cm]")
pl.ylabel("u [m/s]")
pl.grid("on")
pl.tight_layout()
pl.show()

# -------------------------------------
# Synthetic measures
omega = np.matmul(np.transpose(rr), u_gt)
omega = np.matmul(rr, omega)
pl.figure(num=None, figsize=(6, 4))
pl.plot(mesh, omega, label="Noise-free measures")
pl.plot(mesh, omegaF2, label="$\\omega$ ")
pl.legend(loc='best')
pl.xlabel("x [cm]")
pl.ylabel("u [m/s]")
pl.grid("on")
pl.tight_layout()
pl.show()

# -------------------------------------
# Compute reconstruction
u_star_synt = assimilate(mesh, rr, U, omega, nbModes) # synthetic
u_tilde = assimilate(mesh, rr, U, omegaF2, nbModes) # Field2
pl.figure(num=None, figsize=(6, 4))
pl.plot(mesh, u_gt, label="gt", linewidth=1.0)
pl.plot(mesh, u_star_synt, "--", label="$u*(\\omega)$ synt")
pl.plot(mesh, u_tilde, "--", label="$u*(\\omega_{F2})$")
pl.legend(loc='best')
pl.xlabel("x [cm]")
pl.ylabel("u [m/s]")
pl.grid("on")
pl.tight_layout()
pl.show()

# -------------------------------------
# Measure u_tilde
omega_tilde = np.matmul(np.transpose(rr), u_tilde)
omega_tilde = np.matmul(rr, omega_tilde)

# Use expected experimental value to compute R(u_tilde)
meshF2 = np.loadtxt("../../../../data/measures/Doppler_F2/measuresF2_position.txt")
R_u_tilde_coord = data_mean(meshF2, F2, rr, u_tilde)
R_u_tilde = np.zeros(len(mesh))
for i in range(nbMeasures):
  R_u_tilde = R_u_tilde_coord[i]*rr_nn[:,i] + R_u_tilde
omega_unbiased = 2*omega_tilde - R_u_tilde

pl.figure(num=None, figsize=(6, 4))
pl.plot(mesh, u_gt)
pl.plot(mesh, u_tilde)
pl.plot(mesh, omega_tilde, label="$\\Pi_W u_0^*$")
pl.plot(mesh, omegaF2, label="$\\omega$ Field 2")
pl.plot(mesh, R_u_tilde, label="$\\mathbb{E}[R(u_0^*)]$")
pl.plot(mesh, omega_unbiased, label="$\\omega_{unbiased}$")
pl.xlabel("x [cm]")
pl.ylabel("u [m/s]")
pl.legend(loc='best')
pl.grid("on")
pl.tight_layout()
pl.show()

# -------------------------------------
# Compute corrected reconstruction 
u_star = assimilate(mesh, rr, U, omega_unbiased, nbModes) # synthetic
pl.figure(num=None, figsize=(6, 4))
pl.plot(mesh, u_gt, label="$u_{gt}$", linewidth=1.5)
pl.plot(mesh, u_star, linewidth=1.5, label="bPBDW")
#pl.plot(mesh, u_star_synt, "--", label="$u*(\\omega)$ synthetic")
pl.plot(mesh, u_tilde, linewidth=1.0, label="PBDW")
pl.xlabel("x [cm]")
pl.ylabel("u [m/s]")
pl.tight_layout()
pl.legend(loc='best')
pl.grid("on")
pl.show()
