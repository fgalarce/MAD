import numpy as np
import matplotlib.pyplot as pl
from matplotlib import rc
from numpy import linalg as la
import scienceplots
import os

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

firstTime = 0
timeSteps = 10000
nbVertices = 501
dt=1e-3

pl.style.use(['science', 'nature'])
dirResults = os.environ['MAD_RESULTS'] + "/hammer/results/"

density=2000.0
bingham=0.5
viscosity=0.08
length=200
heigth=100.0
diameter=0.1
wave_velocity=1000

Pr=heigth*density*9.81
Ur=1/32*diameter*diameter*Pr/length/viscosity
ca=length/wave_velocity
v = np.loadtxt(dirResults + 'ctrlU.txt');
pl.plot(np.arange(timeSteps)*dt*ca, v[:,0]*Ur, "--", label="$x=0$")
pl.plot(np.arange(timeSteps)*dt*ca, v[:,1]*Ur, label="$x=L/2$")
pl.xlabel("Time [s]")
pl.ylabel("Velocity [m/s]")
pl.xlim([0, timeSteps*dt*ca])
pl.legend()
#pl.tight_layout()
#pl.savefig("pl_dimensional_u.pdf", format="pdf", bbox_inches="tight")
pl.show()

v = np.loadtxt(dirResults + 'ctrlP.txt');
pl.plot(np.arange(timeSteps)*dt*ca, v[:,1]*Pr, "--", label="$x=L/2$")
pl.plot(np.arange(timeSteps)*dt*ca, v[:,2]*Pr, label="$x=L$")
#pl.tight_layout()
pl.xlabel("Time [s]")
pl.ylabel("Pressure [Pa]")
pl.xlim([0, timeSteps*dt*ca])
pl.legend()
#pl.savefig("pl_dimensional_p.pdf", format="pdf", bbox_inches="tight")
pl.show()
