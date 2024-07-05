from tools import *
from assimilation import *

EPSILON = 1e-6
ancho_fig = 6
alto_fig = 4
shut_up = 0

# Parameters
nbSnapshots = 256
nbModes = 26
nbMeasures = 6
pixel_size = 0.1
im_start = np.pi/3
im_end = im_start + pixel_size*nbMeasures

alpha = 0.1
    
# Mesh
mesh = np.arange(0,2.5*np.pi, 0.01)

nbTestCasesA = 8
nbTestCasesT = 8
nbTestCases = nbTestCasesT*nbTestCasesA
error = np.zeros((nbModes-1, nbTestCases))
errorPBDW = np.zeros((nbModes-1, nbTestCases))
iterator = 0

# -------------------------------------
# Manifold generation
A_start = 1.0
A_end = 2.0
T_start = 1.0*np.pi
T_end = 2.0*np.pi
T = np.arange(T_start, T_end, (T_end - T_start)/nbSnapshots)
A = np.arange(A_start, A_end, (A_end - A_start)/nbSnapshots)
#    print(np.shape(T))
manifold = np.zeros((len(mesh), nbSnapshots))
for i in range(nbSnapshots):
  manifold[:,i] = A[i]*np.sin(2*np.pi/T[i]*mesh)
if (shut_up == 0):
  pl.figure(num=None, figsize=(ancho_fig, alto_fig))
  for i in range(nbSnapshots):
    pl.plot(mesh, manifold[:,i])
  pl.xlabel("x")
  pl.tight_layout()
  pl.grid("on")
  pl.xlim([0, 2*np.pi])
#  pl.savefig("2023-07-30_ex1_manifold.pdf", format="pdf", bbox_inches="tight")
  pl.show()
  
# -------------------------------------
# Model reduction
U, S, Vt = np.linalg.svd(manifold, full_matrices=True, compute_uv=True)
mesh = np.arange(0,2*np.pi, 0.01)
U = U[0:len(mesh),:]
S = S[0:len(mesh)]
Vt = Vt[0:len(mesh),:]
if (shut_up == 0):
  pl.figure(num=None, figsize=(ancho_fig, alto_fig))
  for i in range(nbModes):
    pl.plot(mesh, U[:,i], label='mode ' + str(i+1))
  pl.xlim([0, 2*np.pi])
  pl.ylim([-0.1, 0.1])
  pl.tight_layout()
  pl.xlabel("x")
  pl.grid("on")
#  pl.savefig("2023-07-30_ex1_pod.pdf", format="pdf", bbox_inches="tight")
  pl.show()
#      print("is POD orthogonal: ", is_ortho(U))
#      print("is POD normalized: ", is_normal(U))


errorPOD = np.zeros(len(S))
for i in range(len(S)-1):
  errorPOD[i] = np.sqrt(sum(S[i+1:len(S)])) / np.sqrt(sum(S))

# plot model quality
if (shut_up == 0):
  pl.figure(num=None, figsize=(ancho_fig, alto_fig))
  pl.plot(range(1,len(S)+1), errorPOD)
  pl.yscale("log")
  pl.tight_layout()
  pl.xlabel("n")
  pl.grid("on")
#  pl.savefig("2023-07-30_ex1_error_pod.pdf", format="pdf", bbox_inches="tight")
  pl.show()

# -------------------------------------
# Riesz representers
rr = np.zeros((len(mesh), nbMeasures))
if (shut_up == 0):
  pl.figure(num=None, figsize=(ancho_fig, alto_fig))
for i in range(2*nbMeasures):
  if (i % 2 == 0):
    idRiesz = int(i/2)
    pointInside = 0
    for j in range(len(mesh)):
      if (mesh[j] >= im_start + i * pixel_size and mesh[j] < im_start + (i+1) * pixel_size):
        rr[j, idRiesz] = 1.0 
#        print(idRiesz)
    rr[:,idRiesz] = rr[:,idRiesz] / np.linalg.norm(rr[:,idRiesz])
    if (shut_up == 0):
      pl.plot(mesh, rr[:,idRiesz])

for A_gt in np.random.uniform(1,2,nbTestCasesA):
  for T_gt in np.random.uniform(T_start,T_end, nbTestCasesT):
    
#    A_gt = 1.5 
#    T_gt = 1.5*np.pi
    std_dev = A_gt / 100 * 1
          
    # Ground truth
    u_gt = A_gt * np.sin(2*np.pi/T_gt*mesh)
    def u_gt_func(t):
      return A_gt * np.sin(2*np.pi/T_gt*t);
    
    if (shut_up == 0):
      pl.tight_layout()
      pl.grid("on")
      pl.xlabel("x")
      pl.xlim([0, 2*np.pi])
#      pl.savefig("2023-07-30_ex1_rr.pdf", format="pdf", bbox_inches="tight")
      pl.show()
    
#    print("are rr orthogonal: ", is_ortho(rr))
#    print("are rr normalized: ", is_normal(rr))
    
    # -------------------------------------
    # synthetic measures
    omega = np.matmul(np.transpose(rr), u_gt)
    for i in range(nbMeasures):
      omega[i] = omega[i] + np.random.normal(-(omega[i])*alpha, std_dev, 1)
    omega = np.matmul(rr, omega)
    if (shut_up == 0):
      pl.figure(num=None, figsize=(ancho_fig, alto_fig))
      pl.plot(mesh, u_gt)
      pl.plot(mesh, omega)
      pl.tight_layout()
      pl.grid("on")
      pl.xlabel("x")
      pl.xlim([0, 2*np.pi])
#      pl.savefig("2023-07-30_ex1_measures_and_gt.pdf", format="pdf", bbox_inches="tight")
      pl.show()
    
    # -------------------------------------
    # compute a-priori error bound
    G = np.matmul(np.transpose(rr), U[:,0:nbSnapshots])
    # inf-suf constant
    LHS = np.matmul(np.transpose(G), G)
    UU, SS, VVt = np.linalg.svd(LHS, full_matrices=True, compute_uv=True)
    if (shut_up == 0):
      pl.figure(num=None, figsize=(ancho_fig, alto_fig))
      pl.plot(range(1,nbMeasures), SS[0:nbMeasures-1])
      pl.yscale("log")
      pl.tight_layout()
      pl.grid("on")
      pl.xlabel("n")
#      pl.savefig("2023-07-30_ex1_decay.pdf", format="pdf", bbox_inches="tight")
      pl.show()
      
      pl.figure(num=None, figsize=(ancho_fig, alto_fig))
      pl.plot(range(1,nbMeasures-3), errorPOD[0:nbMeasures-4]/SS[0:nbMeasures-4])
      pl.tight_layout()
      pl.yscale("log")
      pl.grid("on")
      pl.xlabel("x")
      pl.xlabel("n")
#      pl.savefig("2023-07-30_ex1_error_pod.pdf", format="pdf", bbox_inches="tight")
      pl.show()
    
    # -------------------------------------
    # normal equations
    #pl.figure(num=None, figsize=(ancho_fig, alto_fig))
    u_star = np.zeros((len(mesh), nbModes))
    error_l2 = np.zeros(nbModes-1)
    for j in range(1,nbModes):
      G = np.matmul(np.transpose(rr), U[:,0:j])
      LHS = np.matmul(np.transpose(G), G)
      RHS = np.matmul(np.transpose(U[:,0:j]), omega)
      c = np.linalg.solve(LHS, RHS)
      for i in range(j):
        u_star[:,j] = u_star[:,j] + U[:,i] * c[i]
      error_l2[j-1] = np.linalg.norm(u_star[:,j] - u_gt) / np.linalg.norm(u_gt)
    errorPBDW[:,iterator] = error_l2
    
    if (shut_up == 0):
      pl.figure(num=None, figsize=(ancho_fig, alto_fig))
      pl.plot(range(1,nbModes), error_l2)
      pl.legend(loc='best')
      pl.yscale("log")
      pl.tight_layout()
      pl.grid("on")
      pl.xlabel("n")
#      pl.savefig("2023-07-30_ex1_rec_l2Error.pdf", format="pdf", bbox_inches="tight")
      pl.show()
    
    # compute a-priori error bound
    G = np.matmul(np.transpose(rr), U[:,0:nbSnapshots])
    LHS = np.matmul(np.transpose(G), G)
    UU, SS, VVt = np.linalg.svd(LHS, full_matrices=True, compute_uv=True)
    if (shut_up == 0):
      pl.figure(num=None, figsize=(ancho_fig, alto_fig))
      # inf-suf constant
      pl.plot(range(1,nbModes), SS[0:nbModes-1])
      pl.xlabel("$n$", fontsize=fontSize)
      pl.ylabel("$\\beta$", fontsize=fontSize)
      pl.yscale("log")
      pl.show()
      
      pl.figure(num=None, figsize=(ancho_fig, alto_fig))
      pl.plot(range(1,nbModes+1), errorPOD[0:nbModes]/SS[0:nbModes])
      pl.xlabel("$n$", fontsize=fontSize)
      pl.ylabel("$\\beta(V_n,W_m)^{-1} \\epsilon_n$", fontsize=fontSize)
      pl.yscale("log")
      pl.show()
    
    # -------------------------------------
    # normal equations
    if (shut_up == 0):
      pl.figure(num=None, figsize=(ancho_fig, alto_fig))
    u_star_uPBDW = np.zeros((len(mesh), nbModes))
    error_l2 = np.zeros(nbModes-1)
    for j in range(1,nbModes):
      eta = np.matmul(np.transpose(rr), u_star[:,j])
      for i in range(nbMeasures):
        position = i*pixel_size/2.0 + im_start
        eta[i] = eta[i] + np.random.normal(eta[i]*alpha, 0.0, 1)
      eta = np.matmul(rr, eta)
    
      G = np.matmul(np.transpose(rr), U[:,0:j])
      LHS = np.matmul(np.transpose(G), G)
      RHS = np.matmul(np.transpose(U[:,0:j]), eta)
      c = np.linalg.solve(LHS, RHS)
      for i in range(j):
        u_star_uPBDW[:,j] = u_star_uPBDW[:,j] + U[:,i] * c[i]
      error_l2[j-1] = np.linalg.norm(u_star_uPBDW[:,j] - u_gt[:]) / np.linalg.norm(u_gt[:])
      if (shut_up == 0):
        pl.plot(mesh, u_star_uPBDW[:,j], label='$u^*$ (' + str(j) + ' modes)')
    
    if (shut_up == 0):
      pl.plot(mesh, u_gt, '--', label='ground truth', linewidth = 2.0)
      pl.legend(loc='best')
      pl.tight_layout()
      pl.grid("on")
      pl.xlabel("x")
      pl.xlim([0, 2*np.pi])
#      pl.savefig("2023-07-30_ex1_reconstruction.pdf", format="pdf", bbox_inches="tight")
      pl.show()
      
      pl.figure(num=None, figsize=(ancho_fig, alto_fig))
      pl.plot(mesh, u_gt, "--", label='ground truth')
      pl.plot(mesh, u_star[:,np.argmin(error_l2)], label='PBDW (' + str(np.argmin(error_l2)+1) + ' modes)')
      pl.plot(mesh, u_star_uPBDW[:,np.argmin(error_l2)], label='bPBDW (' + str(np.argmin(error_l2)+1) + ' modes)')
      pl.legend(loc='best')
      pl.tight_layout()
      pl.grid("on")
      pl.xlabel("x")
      pl.xlim([0, 2*np.pi])
#      pl.savefig("2023-07-30_ex1_reconstruction_comparison.pdf", format="pdf", bbox_inches="tight")
      pl.show()
    error[:,iterator] = error_l2
    iterator = iterator + 1

error_avg = np.zeros(nbModes-1)
error_avgPBDW = np.zeros(nbModes-1)
pl.figure(num=None, figsize=(6, 4.5), dpi=200)
for i in range(nbTestCases):
  error_avg = error_avg + error[:,i]
  error_avgPBDW = error_avgPBDW + errorPBDW[:,i]
  pl.plot(range(1,nbModes), error[:,i], linewidth = 1.0, alpha = 0.2)
  pl.legend(loc='best')
  pl.tight_layout()
  pl.grid("on")
  pl.yscale("log")
  pl.xlabel("n")
error_avg = error_avg * float(1/float(nbTestCases))
pl.plot(range(1,nbModes), error_avg, "k", linewidth = 1.0)
pl.show()

pl.figure(num=None, figsize=(ancho_fig, alto_fig))
#pl.figure(num=None, figsize=(6, 4), dpi=200)
error_avgPBDW = error_avgPBDW * float(1/float(nbTestCases))
np.savetxt("error_avg_" + str(alpha) + ".txt", error_avg)
np.savetxt("error_avgPBDW_" + str(alpha) + ".txt", error_avg)
pl.plot(range(1,nbModes), error_avg, linewidth = 1.0, label="bPBDW" )
pl.plot(range(1,nbModes), error_avgPBDW, linewidth = 1.0, label="PBDW")
print(np.min(error_avg))
pl.legend(loc='best')
pl.tight_layout()
pl.grid("on")
pl.yscale("log")
pl.xlabel("n")
pl.xlabel("n")
#pl.xlim([1, 8])
#pl.savefig("2023-07-30_ex1_reconstruction_error_avg_comp.pdf", format="pdf", bbox_inches="tight")
pl.show()
