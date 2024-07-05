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

fontSize=14
pl.figure(num=None, figsize=(6, 4))
pl.rcParams.update({'font.size': fontSize})
pl.tight_layout()
pl.grid("on")

a = np.loadtxt('./new/p1.csv', skiprows=1, delimiter=",");
b = np.loadtxt('./new/p10.csv', skiprows=1, delimiter=",");
c = np.loadtxt('./new/p25.csv', skiprows=1, delimiter=",");
d = np.loadtxt('./new/p50.csv', skiprows=1, delimiter=",");
e = np.loadtxt('./new/p75.csv', skiprows=1, delimiter=",");
f = np.loadtxt('./new/p100.csv', skiprows=1, delimiter=",");
g = np.loadtxt('./new/p125.csv', skiprows=1, delimiter=",");

#pl.style.use(['science', 'nature'])
pl.plot(a[:,0], a[:,1], label='$\\omega = 1$ Hz')
pl.plot(b[:,0], b[:,1], label='$\\omega = 10$ Hz')
pl.plot(c[:,0], c[:,1], label='$\\omega = 25$ Hz')
pl.plot(d[:,0], d[:,1], label='$\\omega = 50$ Hz')
#pl.plot(e[:,0], e[:,1], label='$\\omega = 75$ Hz')
#pl.plot(f[:,0], f[:,1], label='$\\omega = 100$ Hz')
#pl.plot(g[:,0], g[:,1], label='$\\omega = 125$ Hz')
pl.xlabel("Position [cm]")
pl.ylabel("Pressure [dyn/cm$^2$]")
#pl.yscale("log")
pl.legend(loc='best')
pl.xlim([0,1.15])
pl.axvline(x = 0.333, color = 'k')
pl.axvline(x = 0.666, color = 'k')
#pl.legend(loc='center left', bbox_to_anchor=(1, 0.5))
pl.tight_layout()
#pl.grid("on")
pl.show()

