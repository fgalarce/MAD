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

p = loadtxt('./results/svGramian.txt', skiprows=0);
n=20

plot(range(1,n+1), p, linewidth=2.0)
##legend(loc='best')
grid()
#xticks(fontsize=tickSize, rotation = 0)
#yticks(fontsize=tickSize, rotation = 45)
#ax = axes()
#ax.set_xlabel("n", fontsize=fontSize)
show()
#ax.set_ylabel("", fontsize=fontSize)
#yscale("log")
##savefig(fileName)
#clf()

#call('open ' + fileName, shell=True)
