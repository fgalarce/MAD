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

fontSize=25
tickSize=20
figure(num=None, figsize=(9, 6))
rcParams.update({'font.size': fontSize})

p = loadtxt('./EVD.txt', skiprows=0);
fileName = 'EV_mds.pdf'

plot(p[0,:]/p[0,0], linewidth=2.0)
legend(loc='best')
grid()
xticks(fontsize=tickSize, rotation = 0)
yticks(fontsize=tickSize, rotation = 45)
ax = axes()
ax.set_xlabel("$n$", fontsize=fontSize)
ax.set_ylabel("$\lambda$", fontsize=fontSize)
yscale("log")
show()
#savefig(fileName)
#clf()
#
#call('open ' + fileName, shell=True)
