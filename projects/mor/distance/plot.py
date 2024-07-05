from numpy import *
from matplotlib import *
from matplotlib.pyplot import *
from subprocess import *
from matplotlib import rc

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

rc('text', usetex=True)
style.use('ggplot')

figure(num=None, figsize=(10, 6))
rcParams.update({'font.size': 16})

nbMeshes=25
dA = zeros(nbMeshes)
dBC = zeros(nbMeshes)
dC = zeros(nbMeshes)
dFS = zeros(nbMeshes)
dM = zeros(nbMeshes)
dP = zeros(nbMeshes)
dProj = zeros(nbMeshes)
dS = zeros(nbMeshes)
dG = zeros(nbMeshes)
dO = zeros(nbMeshes)
dH = zeros(nbMeshes)

for i in range(1,nbMeshes+1):
  p = loadtxt('./results/distance.' + str(i) + ".txt", skiprows=0);
  q = loadtxt('./results/distance_olga' + str(i) + ".txt", skiprows=0);
  dO[i-1] = q;
  dA[i-1] = p[0];
  dBC[i-1] = p[1]
  dC[i-1] = p[2]
  dFS[i-1] = p[3]
  dM[i-1] = p[4]
  dP[i-1] = p[5]
  dProj[i-1] = p[6]
  dS[i-1] = p[7]
  dG[i-1] = p[8]
  dH[i-1] = p[9]

print(dH)
fileName = 'distance.pdf'
isort = argsort(dH)
#plot(range(1,nbMeshes+1), dA[isort]/max(dA), label='Asimov', linewidth=2.0)
#plot(range(1,nbMeshes+1), dO[isort]/max(dO), label='Olga', linewidth=2.0)
plot(range(1,nbMeshes+1), dH[isort], label='Hausdorff', linewidth=2.0)
#plot(range(1,nbMeshes+1), dBC/max(dBC), label='Binet-Cauchy', linewidth=2.0)
#plot(range(1,nbMeshes+1), dC[isort]/max(dC), label='Chordal', linewidth=2.0)
#plot(range(1,nbMeshes+1), dFS/max(dFS), label='Fubini-Study', linewidth=2.0)
#plot(range(1,nbMeshes+1), dM[isort]/max(dM), label='Martin', linewidth=2.0)
#plot(range(1,nbMeshes+1), dP/max(dP), label='Procrustes', linewidth=2.0)
#plot(range(1,nbMeshes+1), dProj[isort]/max(dProj), label='Projection', linewidth=2.0)
#plot(range(1,nbMeshes+1), dS[isort]/max(dS), label='Spectral', linewidth=2.0, linestyle='--')
#plot(range(1,nbMeshes+1), dG[isort]/max(dG), label='Grassmanian', linewidth=2.0)

legend(loc='best')
grid()
xticks(fontsize=16, rotation = 0)
yticks(fontsize=16, rotation = 45)
ax = axes()
ax.set_xlabel("i")
ax.set_ylabel("Relative distance")
savefig(fileName)
clf()

call('open ' + fileName, shell=True)
