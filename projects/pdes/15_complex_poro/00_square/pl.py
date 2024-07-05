import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc
import seaborn as sns
sns.set_theme()
fontSize=18
pl.figure(num=None, figsize=(8, 6), dpi=200)
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

p = np.loadtxt('./convergence2.txt', skiprows=0);
pl.plot(p[:,1],p[:,2], "*", color="k")
pl.plot(p[:,1],p[:,2], label="$||\\vec{u} - \\vec{u}^h||$", color="k")
pl.plot(p[:,1],p[:,3], "*", "r")
pl.plot(p[:,1],p[:,3], label="$||u_{re} - u_{re}^h||$", color="r")
pl.plot(p[:,1],p[:,4], "*", color="b")
pl.plot(p[:,1],p[:,4], label="$||u_{im} - u_{im}^h||$", color="b")
pl.plot(p[:,1],p[:,5], "*", color="g")
pl.plot(p[:,1],p[:,5], label="$||p_{re} - p_{re}^h||$", color="g")
pl.plot(p[:,1],p[:,6], "*", color="y")
pl.plot(p[:,1],p[:,6], label="$||p_{im} - p_{im}^h||$", color="y")
pl.plot(p[:,1],p[:,7], "*", color="m")
pl.plot(p[:,1],p[:,7], label="$||\\phi_{re} - \\phi_{re}^h||$", color="m")
pl.plot(p[:,1],p[:,8], "*", color="c")
pl.plot(p[:,1],p[:,8], label="$||\\phi_{im} - \\phi_{im}^h||$", color="c")

pl.legend(loc='best', fontsize=10, bbox_to_anchor=(1, 0.5))
pl.xlabel("h [cms]", fontsize=fontSize)
#pl.ylabel("$|| u - u_h ||_U$", fontsize=fontSize)
pl.yscale("log")
pl.xlim(max(p[:,1]), min(p[:,1] ))
pl.tight_layout()
#pl.grid()
pl.show()
