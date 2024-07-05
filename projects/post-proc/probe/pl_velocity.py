import numpy as np
import fileinput
import matplotlib.pyplot as pl
from matplotlib import rc
import argparse
parser = argparse.ArgumentParser(description='parser')
parser.add_argument('-n', action="store", dest='simNumber', default=1)
args = parser.parse_args()

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

simNumber =int(args.simNumber)

p = np.loadtxt('./results/sim' + wildcard(simNumber) + '/velocity_at_probe.txt', skiprows=0);

dt = 0.01
pl.plot(dt*np.arange(100), p[:,0], label='x')
pl.plot(dt*np.arange(100), p[:,1], label='y')
pl.legend(loc='best')
pl.xlabel("Time [sec]")
pl.ylabel("Velocity [m/s]")
pl.show()
