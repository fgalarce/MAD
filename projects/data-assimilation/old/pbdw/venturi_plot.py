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

p = loadtxt('./results/error_pbdw.txt', skiprows=1);
fileName = 'pbdw_venturi.pdf'

plot(p[:,0]/p[:,1], linewidth=3.0)
legend(loc='best')
grid()
xticks(fontsize=14, rotation = 0)
yticks(fontsize=14, rotation = 45)
ax = axes()
ax.set_xlabel("Time iteration")
ax.set_ylabel("$H^1$ error")
savefig(fileName)
clf()

call('open ' + fileName, shell=True)
