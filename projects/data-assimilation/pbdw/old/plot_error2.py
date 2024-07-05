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

timeStep=.02
nbSimulations=16

error = zeros(5)
error_ps = zeros((nbSimulations, 5))
k = 0
for n in range(10,60,10):
  avgError = 0
  for i in range(nbSimulations):
    p = loadtxt('./results/sim' + wildcard(i) + '/n' + str(n) + '/error.txt', skiprows=0)
    avgError = avgError + sqrt(sum(p[:,0]*p[:,0]))/sqrt(sum(p[:,1]*p[:,1]))
    error_ps[i][k] = sqrt(sum(p[:,0]*p[:,0]))/sqrt(sum(p[:,1]*p[:,1]))
  avgError = avgError/nbSimulations
  error[k] = avgError
  k = k + 1

for i in range(nbSimulations):
  plot(range(10,60,10), error_ps[i], alpha=0.25)

plot(range(10,60,10), error, linewidth=2.5)
#legend(loc='best')
#grid()
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
