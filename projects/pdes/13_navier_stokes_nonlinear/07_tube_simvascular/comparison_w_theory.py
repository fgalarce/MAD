import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc
import scienceplots
import os
dir_results = os.environ['MAD_RESULTS']

# Parameters
Q0 = 10.0
Rd = 1.0
Rp = 0.1
C  = 1.0/(4.0*np.pi)

# Time
dt = 1e-3
t = np.arange(0.0, 500*dt, dt)

# Theoretical pressure
ts = t/(Rd*C)
Pex = Rd*Q0*( (Rp/Rd + 0.5 )*np.sin(ts/2)**2 + 0.25*(1.0 - np.exp(-ts) - np.sin(ts)) )

# Import pressure from simulation
P_sim = np.loadtxt(dir_results + "/ns_tube_simvascular/wkssl_pressure.txt")
P_sim2 = np.loadtxt(dir_results + "/ns_tube_simvascular_old/wkssl_pressure.txt")
t_sim = np.linspace(dt, len(P_sim)*dt, len(P_sim))/(Rd*C)
t_sim2 = np.linspace(dt, len(P_sim2)*dt, len(P_sim2))/(Rd*C)


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


pl.plot(ts/np.pi, Pex, label="Theory")
pl.plot(t_sim/np.pi, P_sim, label = "Simulation")
pl.plot(t_sim2/np.pi, P_sim2, label = "Simulation Old")
pl.xlabel("Time [s]")
pl.ylabel("$p_d + p_p$")
pl.legend(loc='best')
pl.show()
