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

dx = 0.03
x = np.arange(-0.25, 0.25 + dx, dx)
n = 0.5
u_max = (2*n+1)/(n+1) * 0.8
u_an = u_max * (1.0 - (x/0.25)**((n+1.0)/n))**((n+1.0)/n)
pl.plot(x, u_an, label='analytic', linewidth=2.0)
#u = np.loadtxt('numerical.txt', skiprows=1);
#pl.plot(x, u, label='numerical', linewidth=2.0)

#pl.legend(loc='best')
#pl.xlabel("<++>", fontsize=fontSize)
#pl.ylabel("<++>", fontsize=fontSize)
#pl.yscale("log")
#pl.tight_layout()
#pl.grid()
pl.show()
