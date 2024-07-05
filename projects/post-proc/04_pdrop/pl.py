import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc
import scienceplots
import os
dir_results = os.environ['MAD_RESULTS']

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

p = np.loadtxt('13kv_pdrop.txt', skiprows=0)[:,:];
q = np.loadtxt('24kv_pdrop.txt', skiprows=0)[:,:];
r = np.loadtxt('36kv_pdrop.txt', skiprows=0)[:,:];
s = np.loadtxt('57kv_pdrop.txt', skiprows=0)[:,:];
t = np.loadtxt('96kv_pdrop.txt', skiprows=0)[:,:];

pl.style.use(['science', 'nature'])
pl.plot(p[:,0], linewidth=1.0, label="13kv")
pl.plot(q[:,0], linewidth=1.0, label="24kv")
pl.plot(r[:,0], linewidth=1.0, label="36kv")
pl.plot(s[:,0], linewidth=1.0, label="57kv")
pl.plot(t[:,0], linewidth=1.0, label="96kv")
pl.xlabel("Time step")
pl.ylabel("Pressure CGS")
#pl.yscale("log")
pl.legend(loc='best')
#pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))
pl.tight_layout()
#pl.grid("on")
pl.show()

#pl.subplot(221)
#pl.cla()
#pl.subplot(222)
#pl.cla()
#pl.subplot(223)
#pl.cla()
