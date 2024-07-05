from numpy import *
from matplotlib import *
from matplotlib.pyplot import *
from subprocess import *
from matplotlib import rc
rc('text', usetex=True)

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

figure(num=None, figsize=(9, 7), dpi=80, facecolor='w', edgecolor='k')
rcParams.update({'font.size': 20})
latexSize=20
tickSize=20
style.use('ggplot')

bound = 31.659

pmean = zeros((30,))
for i in range(16):
  p = loadtxt('/local/fgalarce/4d-flow-ultrasound/data/reconstruction/up_H1/sim' + wildcard(i) + "/pdrop.txt", skiprows=1);
  pmean = pmean + p[0:30,2]

pmean = pmean/16

time = arange(p[0:30,2].size)*0.004
plot(time, bound/pmean, linewidth=1.5)


fileName = 'pdrop_bound_gamma1.pdf'

#yscale('log')
legend(loc='best')
grid()
xticks(fontsize=tickSize, rotation = 0)
yticks(fontsize=tickSize, rotation = 45)
ax = axes()
ax.set_xlabel("Time [sec]", fontsize=latexSize)
ax.set_ylabel("$e(\delta_p)$", fontsize=latexSize)
savefig(fileName)
clf()

call('evince ' + fileName, shell=True)
