from subprocess import *
from numpy import *
from multiprocessing import Pool
from fileinput import *

def getParameter(dataFilePath, variableName):
  datafile = open(dataFilePath, 'r')
  line = datafile.readline()
  while line:
    if (line.split("=")[0] == variableName):
      return line.split("=")[1]
    line = datafile.readline()

def setParameter(parameter, value, fileName):
  datafile = open(fileName, 'r+')
  succes = False
  for line in input( fileName ):
    if (line.split("=")[0] == parameter):
      datafile.write(line.replace(line[0:len(line)-1], line.split("=")[0] + "=" +  str(value)))
      succes = True
    else:
      datafile.write(line)
  datafile.close()
  if not succes:
    exit("ERROR: Parameter " + parameter + " not found in " + fileName)

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

# Define the parser
import argparse
parser = argparse.ArgumentParser(description='parser')
parser.add_argument('-n', action="store", dest='n_procs', default=1)
args = parser.parse_args()

n_procs =int(args.n_procs)

nbMeshes = 64
h_voxel = 0.09
noiseIt = 100
snr = 0.3

def loop(j):
  for i in range(noiseIt):
    g = loadtxt('../../data/venturi/mesh/v' + wildcard(j) + "/g.txt")
    # Noise as a fraction of voxel volume
    g_noise = g + random.normal(0, (h_voxel*h_voxel*h_voxel)*snr, 1800) 
    savetxt('../../data/venturi/mesh/v' + wildcard(j) + '/g' + wildcard(i) + '.txt', g_noise) 
    
if n_procs == 1:
  for j in range(nbMeshes):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(nbMeshes)) 
