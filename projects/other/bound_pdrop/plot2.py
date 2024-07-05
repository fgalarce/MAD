from numpy import *
from matplotlib import *
from matplotlib.pyplot import *
from subprocess import *

mmHg = 760.0 / (10.0 * 101325.0)
mu = 23.7617

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

from matplotlib import rc
rc('text', usetex=True)

figure(num=None, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
rcParams.update({'font.size': 30})
latexSize=30
tickSize=30
style.use('ggplot')

plot_type = 1

eps_n = zeros(250)
max_eps_n = zeros(250)
i=5

for j in range(i):
  for k in range(i):
    q = loadtxt("./" + str(i) + "_" + str(i) + "/" + str(j) + "_" + str(k) + "/window_time" + str(j) + "_HR" +str(k) + "/eigenValues.txt" )

    for n in range(250):
      eps_n[n] = sqrt( sum(q[n:q.size]) / sum(q))
      if (eps_n[n] > max_eps_n[n]):
        max_eps_n[n] = eps_n[n]

n = 40
p = zeros(n-1)
for n in range(1,40):
  p[n-1] = max(loadtxt('./results/mu' + str(n) + '.txt', skiprows=0));
plot(max_eps_n[10:39]*p[10:39], linewidth=1.5)

yscale('log')
grid()
xticks(fontsize=tickSize, rotation = 0)
yticks(fontsize=tickSize, rotation = 45)
#ax = axes()

fileName = './eps_n_kappa.pdf'
savefig(fileName)
clf()
