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
nbSimulations=16

for i in range(nbSimulations):
  p = loadtxt('./results/sim' + wildcard(i) + '/n20/error.txt', skiprows=0)
  timeStep = 0.02
  nbIterations = len(p[:,0])
  
  # H1 error
  plot(arange(nbIterations)*timeStep, p[:,0]/sqrt(sum(p[:,1]*p[:,1])), linewidth=1.0)
  # gradient error
#  plot(arange(nbIterations)*timeStep, p[:,3]/sqrt(sum(p[:,2]*p[:,2])), linewidth=1.0)
  # L2 error
#  plot(arange(nbIterations)*timeStep, (p[:,0] - p[:,3])/sqrt(sum(-p[:,2]*p[:,2] + p[:,1]*p[:,1] )), linewidth=1.0)
  
#legend(loc='best')
#grid()
##xticks(fontsize=tickSize, rotation = 0)
##yticks(fontsize=tickSize, rotation = 45)
#ax = axes()
#ax.set_xlabel("time", fontsize=fontSize)
yscale("log")
show()
#ax.set_ylabel("", fontsize=fontSize)
##savefig(fileName)
#clf()

#call('open ' + fileName, shell=True)
