import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc
import scienceplots

pl.style.use(['science', 'nature'])
#pl.style.use(['science'])

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

timeSteps=10000
dt=1e-3*L/a
time= dt*np.arange(timeSteps)

q = np.loadtxt('./results_1000/ctrlP.txt', skiprows=0)[0:timeSteps];
r = np.loadtxt('./results_500/ctrlP.txt', skiprows=0)[0:timeSteps];
s = np.loadtxt('./results_200/ctrlP.txt', skiprows=0)[0:timeSteps];
t = np.loadtxt('./results_100/ctrlP.txt', skiprows=0)[0:timeSteps];

pl.plot(time, q[:,2]*pr, linewidth=2.0, label='$\Delta x^*= 1 \\times 10^{-3}$')
pl.plot(time, r[:,2]*pr, linewidth=1.5, label='$\Delta x^*= 2 \\times 10^{-3}$')
pl.plot(time, s[:,2]*pr, "k", label='$\Delta x^* = \\times 5 \\times 10^{-3}$')
pl.plot(time, t[:,2]*pr, "--", label='$\Delta x^* = \\times 1 \\times 10^{-2}$')
#pl.xlim([-0.01,0.02])
#pl.ylim([1.2*10e6,2.7*10e6])
pl.xlabel("Time [s]")
pl.ylabel("Pressure [Pa]")
pl.legend(loc="upper right")
#pl.tight_layout()
pl.savefig("convergence_space.pdf", format="pdf", bbox_inches="tight")
pl.show()
