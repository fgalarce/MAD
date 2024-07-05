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

fontSize=35
tickSize=25
figure(num=None, figsize=(9, 6))
rcParams.update({'font.size': fontSize})

nbModes = 20
p = zeros(nbModes-3)
q = zeros(nbModes-3)
r = zeros(nbModes-3)

for n in range(3,nbModes):
  p[n-3] = loadtxt('./results_under_almost_Hausdorff/haus' + str(n) + ".txt", skiprows=0);
  q[n-3] = loadtxt('./results_super_almost_Hausdorff/haus' + str(n) + ".txt", skiprows=0);
  r[n-3] = loadtxt('./results_almost_Hausdorff/haus' + str(n) + ".txt", skiprows=0);

fileName = 'super_sampling_almost_Haus.pdf'

plot(range(3,nbModes),p, label='800 (32)', linewidth=2.0)
plot(range(3,nbModes),r, label='1600 (64)', linewidth=2.0)
plot(range(3,nbModes),q, label='12800 (512)', linewidth=2.0)
grid()
xticks(fontsize=tickSize, rotation = 0)
yticks(fontsize=tickSize, rotation = 45)
ax = axes()
ax.set_xlabel("$n$", fontsize=fontSize)
ax.set_ylabel("$d(V_n^i, V_n^j)$", fontsize=fontSize)
legend(loc='best', fontsize=tickSize)
show()
