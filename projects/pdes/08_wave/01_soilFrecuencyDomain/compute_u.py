import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc
import seaborn as sns
sns.set_theme()
fontSize=18
pl.figure(num=None, figsize=(6, 4), dpi=200)
pl.rcParams.update({'font.family': 'serif'})
pl.rcParams.update({'font.size': fontSize})

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

a = np.loadtxt('armenia_1999.txt');

nbIterations = len(a[:,0])
d0 = 0.0
v0 = 0.0
dt=5e-3
d = np.zeros(nbIterations)
for i in range(nbIterations):
  d1 = d0 + dt*(v0 + dt*a[i,1])
  d[i] = d1
  v0 = (d1 - d0)/dt

#pl.plot(a[:,0], a[:,1], linewidth=0.5)
pl.plot(a[:,0], d, linewidth=0.5)
np.savetxt("./disp_armenia.txt", d)

#pl.legend(loc='best')
pl.xlabel("t", fontsize=fontSize)
pl.ylabel("u", fontsize=fontSize)
#pl.yscale("log")
#pl.tight_layout()
#pl.grid()
pl.show()
