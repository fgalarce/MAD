from numpy import *
from matplotlib import *
from matplotlib.pyplot import *
from subprocess import *

from matplotlib import rc
rc('text', usetex=True)

figure(num=None, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
rcParams.update({'font.size': 20})
latexSize=20
tickSize=20
style.use('ggplot')

nbMeshes=60

a = zeros((nbMeshes,nbMeshes))

for i in range(nbMeshes):
  for j in range(nbMeshes):
    a[i][j] = loadtxt('./training/e' + str(i) + "_" + str(j), skiprows=0);

savetxt("E.txt", a)
