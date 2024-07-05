import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc
import scienceplots
import os

pl.style.use(['science', 'nature'])

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

L=200
a=1000
rho = 2000
g = 9.81
h = 100
pr = rho * g * h
factor=2
timeSteps=factor*500


dir_results = os.environ['MAD_RESULTS'] + "/hammer/plot_paper/"

dt=1e-2*L/a
p = np.loadtxt(dir_results + '/results2/ctrlP.txt', skiprows=0)[0:timeSteps];
timep= dt*np.arange(timeSteps)

timeSteps=factor*5000
dt=1e-3*L/a
q = np.loadtxt(dir_results + '/results3/ctrlP.txt', skiprows=0)[0:timeSteps];
timeq= dt*np.arange(timeSteps)

timeSteps=factor*50000
dt=1e-4*L/a
r = np.loadtxt(dir_results + '/results4/ctrlP.txt', skiprows=0)[0:timeSteps];
timer= dt*np.arange(timeSteps)

timeSteps=factor*500000
dt=1e-5*L/a
s = np.loadtxt(dir_results + '/results5/ctrlP.txt', skiprows=0)[0:timeSteps];
times = dt*np.arange(timeSteps)

pl.plot(timep, p[:,1]*pr, label='$\Delta t=2 \\times 10^{-3}$ s')
pl.plot(timeq, q[:,1]*pr, label='$\Delta t=2 \\times 10^{-4}$ s')
pl.plot(timer, r[:,1]*pr, "--", color="black", label='$\Delta t=2 \\times 10^{-5}$ s', linewidth=2.0)
pl.plot(times, s[:,1]*pr, label='$\Delta t=2 \\times 10^{-6}$ s')
pl.xlabel("Time [s]")
pl.ylabel("Pressure [Pa]")
pl.legend()
pl.tight_layout()
pl.savefig("convergence_time.pdf", format="pdf", bbox_inches="tight")
pl.show()
