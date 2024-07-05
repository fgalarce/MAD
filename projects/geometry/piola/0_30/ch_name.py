from subprocess import *
from numpy import *
from multiprocessing import Pool
import fileinput
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

for i in range(20):

  call("mv Pv." + wildcard(i) + ".vct velocity_mode." + wildcard(i) + ".vct", shell=True)
  
