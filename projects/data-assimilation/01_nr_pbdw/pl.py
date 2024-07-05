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

p = np.loadtxt('./error_avg_0.txt', skiprows=0);
pl.plot(range(1,26), p, label='$\\alpha=0$', linewidth=1.0)

p = np.loadtxt('./error_avg_0.05.txt', skiprows=0);
pl.plot(range(1,26), p, label='$\\alpha=0.05$', linewidth=1.0)

p = np.loadtxt('./error_avg_0.1.txt', skiprows=0);
pl.plot(range(1,26), p, label='$\\alpha=0.1$', linewidth=1.0)

p = np.loadtxt('./error_avg_0.2.txt', skiprows=0);
pl.plot(range(1,26), p, label='$\\alpha=0.2$', linewidth=1.0)

pl.legend(loc='best')
pl.xlabel("$n$", fontsize=fontSize)
pl.ylabel("$e(n)$", fontsize=fontSize)
pl.yscale("log")
pl.tight_layout()
pl.grid("on")
pl.savefig("2023-07-30_ex1_error_vs_alpha.pdf", format="pdf", bbox_inches="tight")
pl.show()
