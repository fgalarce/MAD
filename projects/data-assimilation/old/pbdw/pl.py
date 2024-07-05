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
#style.use('fivethirtyeight')

fontSize=25
tickSize=20
figure(num=None, figsize=(9, 6))
rcParams.update({'font.size': fontSize})

fileName = 'pbdw.pdf'

p = loadtxt('./venturi_piola/error_pbdw.txt', skiprows=1);
time = arange(0,0.5, 0.02)
#plot(time, sqrt(p[:,0]*p[:,0] - p[:,3]*p[:,3])/(sqrt(p[:,1]*p[:,1] - p[:,4]*p[:,4])),  linewidth=2.0)
plot(time,p[:,3]/p[:,4], linewidth=2.0)
#legend(loc='best')
grid()
xticks(fontsize=tickSize, rotation = 0)
yticks(fontsize=tickSize, rotation = 45)
ax = axes()
ax.set_xlabel("Time [sec]", fontsize=fontSize)
ax.set_ylabel("Relative L2 error in gradient", fontsize=fontSize)
show()
#savefig(fileName)
#clf()

call('open ' + fileName, shell=True)


