import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc
import os
import seaborn as sns
sns.set_theme()
fontSize=18
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

p = np.loadtxt('./error_time.txt', skiprows=5)
p=p[50:len(p)];
r = np.loadtxt('./error_time_pbdw.txt', skiprows=5);
r=r[50:len(r)];
q = np.ones(len(p))*0.027684350014556398
time=np.arange(0, 7.4*len(p), 7.4)

pl.figure(num=None, figsize=(4, 3))
pl.plot(time, p, label="bPBDW")
pl.plot(time, r, label="PBDW")
#pl.plot(q)
pl.xlabel("Time [ms]")
pl.ylabel("$|| u - u_{\\text{ref}} || / || u_{\\text{ref}} ||$")
pl.tight_layout()
pl.grid("on")
pl.legend(loc='best', fontsize=10)
pl.show()

print(1.0/len(p)*np.sum(r))
print(np.max(r))
print(np.min(r))
print(np.std(r))
