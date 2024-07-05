import sys
sys.path.insert(1, '../../../io/')
sys.path.insert(1, './')
from tools import *
import os, argparse, sys, subprocess, imageio, time, medit, inputOutput
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from torch.autograd import grad
from pykeops.torch.kernel_product.formula import *

# Torch type and device
use_cuda = torch.cuda.is_available()
torchdeviceId = torch.device('cuda:0') if use_cuda else 'cpu'
torchdtype = torch.float32

# PyKeOps counterpart
KeOpsdeviceId = torchdeviceId.index  # id of Gpu device (in case Gpu is  used)
KeOpsdtype = torchdtype.__str__().split('.')[1]  # 'float32'

parser = argparse.ArgumentParser(description='parser')
parser.add_argument('-source', action="store", dest='fileSource', default='../../data/mesh/caro/caro_sten_5kv.mesh')
parser.add_argument('-dimension', action="store", dest='dimension', default=3)
parser.add_argument('-target', action="store", dest='fileTarget', default='../../data/mesh/caro/caro_5kv.mesh')
parser.add_argument('-scale', type=float, action="store", dest='scaleFactor', default=0.1)
parser.add_argument('-correct_fpc', type=int, action="store", dest='correct_fpc', default=0)
parser.add_argument('-transversal_size', type=float, action="store", dest='transversal_size', default=0.1)
parser.add_argument('-name', action="store", dest='name', default='lddmm')
parser.add_argument('-save_folder', action="store", dest='save_folder', default='./results/')
args = parser.parse_args()
name = args.name
save_folder = args.save_folder
dimension = int(args.dimension)
transversal_size = args.transversal_size

# Load the dataset
fileSource, extSource = os.path.splitext(args.fileSource)
fileTarget, extTarget = os.path.splitext(args.fileTarget)
if (extSource != extTarget):
  print('ERROR: source and target meshes should have the same format.')
  sys.exit()
else:
  ext = extSource;

# optimize only on this boundary label
labels = np.array([4])

if (ext == '.geo'):
  VS, FS, VS_3d, FS_3d, v_map_S = inputOutput.importEnsight(fileSource)
  VT, FT, VT_3d, FT_3d, v_map_T = inputOutput.importEnsight(fileTarget)
elif(ext == '.mesh'):
  VS, FS, VS_3d, FS_3d, v_map_S, pointsId_lddmm_s = medit.importMedit(fileSource, args.scaleFactor, dimension, labels)
  VT, FT, VT_3d, FT_3d, v_map_T, pointsId_lddmm_t = medit.importMedit(fileTarget, args.scaleFactor, dimension, labels)
elif(ext == '.pt'):
  VS = torch.load(fileSource + "_vertices.pt")
  FS = torch.load(fileSource + "_triangles.pt")
  VT = torch.load(fileTarget + "_vertices.pt")
  FT = torch.load(fileTarget + "_triangles.pt")

q0 = VS.clone().detach().to(dtype=torchdtype, device=torchdeviceId).requires_grad_(True)
VT = VT.clone().detach().to(dtype=torchdtype, device=torchdeviceId)
FS = FS.clone().detach().to(dtype=torch.long, device=torchdeviceId)
FT = FT.clone().detach().to(dtype=torch.long, device=torchdeviceId)

# - - - - - -  - - LDDMM - - - - - - - - - - - -
from kernel import *
from ode import *

def Shooting(p0, q0, K, nt=10, Integrator=RalstonIntegrator()):
    return Integrator(HamiltonianSystem(K), (p0, q0), nt)

# lddmm loss 
std_dev=0.05
gamma_=1.0
# L-BFGS algorithm 
learning_rate=0.02
max_iterations=20
tol_change=1e-15
tol_grad=1e-9
if (args.correct_fpc == 1):
  max_iterations=int(max_iterations * 0.1 / transversal_size) + 1
  print("max iterations:", max_iterations)

line_function='strong_wolfe' # 'strong_wolfe' or None
nt=1

sigma = torch.tensor([std_dev], dtype=torchdtype, device=torchdeviceId) # this works for the tetra

def correctFPC(lddmm_field, V, transversal_size):
  Center = np.array([0.25, 0.25, 0.0])
  Vnumpy = V.detach().cpu().numpy()
  for i in range(len(lddmm_field)):
    RC = Center - Vnumpy[i]
    corrector = np.linalg.norm(Center[1] - Vnumpy[i][1]) / np.linalg.norm(Center[1] - Center[1] - transversal_size)
    if np.dot(lddmm_field[i], RC) < 0:
      lddmm_field[i] = -lddmm_field[i] * corrector
  return lddmm_field

def LDDMMloss(K, dataloss, gamma=gamma_):
    def loss(p0, q0):
        p,q = Shooting(p0, q0, K)[-1]
        return gamma * Hamiltonian(K)(p0, q0) + dataloss(q)
    return loss

def lossVarifoldSurf(FS, VT, FT, K):
    # VT: vertices coordinates of target surface, 
    # FS,FT : Face connectivity of source and target surfaces
    # K kernel
    def CompCLNn(F, V):
        if dimension == 2:
            V0, V1  = V.index_select(0, F[:, 0]), V.index_select(0, F[:, 1])
            normals = zeros((size(V0,0), 3), dtype=float32);
            for i in range(size(normals,0)):
              normals[i][0] = - V0.detach().cpu().numpy()[i][1] + V1.detach().cpu().numpy()[i][1];
              normals[i][1] = - V1.detach().cpu().numpy()[i][0] + V0.detach().cpu().numpy()[i][0];
              normals[i][2] = 0.0;
              normalization = sqrt(normals[i][0]*normals[i][0] + normals[i][1]*normals[i][1])
              normals[i] = normals[i]/normalization;
            C, N = 0.5 * (V0 + V1), torch.from_numpy(normals)
        else:
            V0, V1, V2 = V.index_select(0, F[:, 0]), V.index_select(0, F[:, 1]), V.index_select(0, F[:, 2])
            C, N = .333 * (V0 + V1 + V2), .5 * torch.cross(V1 - V0, V2 - V0)

        L = (N ** 2).sum(dim=1)[:, None].sqrt()
        return C, L, N / L
    
    CT, LT, NTn = CompCLNn(FT, VT)
    cst = (LT * K(CT, CT, NTn, NTn, LT)).sum()
    
    def loss(VS):
        CS, LS, NSn = CompCLNn(FS, VS)
        geometric_loss = cst + (LS * K(CS, CS, NSn, NSn, LS)).sum() - 2 * (LS * K(CS, CT, NSn, NTn, LT)).sum()
        return geometric_loss
    return loss

# Define data attachment and LDDMM functional
dataloss = lossVarifoldSurf(FS, VT, FT, GaussLinKernel(sigma=sigma))
Kv = GaussKernel(sigma=sigma)
loss = LDDMMloss(Kv, dataloss)

p0 = torch.zeros(q0.shape, dtype=torchdtype, device=torchdeviceId, requires_grad=True)
optimizer = torch.optim.LBFGS([p0], lr=learning_rate, max_iter=max_iterations, tolerance_change=tol_change, tolerance_grad=tol_grad, line_search_fn=line_function)
print('   LDDMM optimzation. pyKeOps + pyTorch')
start = time.time()

def closure():
    optimizer.zero_grad()
    L = loss(p0, q0)
    print('OPTI: Loss', L.detach().cpu().numpy())
    fileLoss = open(save_folder + "/loss_" + name + ".txt", "a+")
    fileLoss.write(str(L.detach().cpu().numpy()) + "\n")
    L.backward() # Compute gradient of L
    return L

nt_shoot = 100

fileLoss = open(save_folder + "/loss_" + name + ".txt", "w")
for i in range(nt):
  print('\nOPTI: it', i)
  optimizer.step(closure)

  listpq = Shooting(p0, q0, Kv, nt_shoot)
  inputOutput.writeCase(save_folder, name, "lddmm", nt)
#  if i == nt - 1:
#    disp_lddmm = ( listpq[nt_shoot-1][1].detach().cpu().numpy() - q0.detach().cpu().numpy() ) / args.scaleFactor
#    inputOutput.writeState(save_folder, name + "." + wildcard(0), VS_3d, disp_lddmm, pointsId_lddmm_s, 0)
  disp_lddmm = ( listpq[nt_shoot-1][1].detach().cpu().numpy() - q0.detach().cpu().numpy() ) / args.scaleFactor
  if args.correct_fpc == 1:
    disp_lddmm = correctFPC(disp_lddmm, VS, transversal_size)
  inputOutput.writeState(save_folder, name + "." + wildcard(i), VS_3d, disp_lddmm, pointsId_lddmm_s, 0)

fileLoss.close()
print(' - - - - - - - - - - - - - - - - - - - - - - - - - - -')
print('Optimization (L-BFGS) time: ', round(time.time() - start, 2), ' seconds')
