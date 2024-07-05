from subprocess import *
from numpy import *
from sys import *
from multiprocessing import Pool
import fileinput
import argparse

parser = argparse.ArgumentParser(description='parser')
parser.add_argument('-n', action="store", dest='n_procs', default=1)
args = parser.parse_args()

n_procs = int(args.n_procs)

def getParameter(dataFilePath, variableName):
  datafile = open(dataFilePath, 'r')
  line = datafile.readline()
  while line:
    if (line.split(" = ")[0] == variableName):
      return line.split("=")[1]
    line = datafile.readline()

def setParameter(parameter, value, fileName):
  datafile = open(fileName, 'r+')
  succes = False
  for line in fileinput.input( fileName ):
    if (line.split("=")[0] == parameter):
      datafile.write(line.replace(line[0:len(line)-1], line.split("=")[0] + "=" +  str(value)))
      succes = True
    else:
      datafile.write(line)
  if not succes:
    datafile.write(parameter +"=" +  str(value))
  datafile.close()
  
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

nbMeshes = 60
#noiseIterations = 100

def loopSim(i):

  dirResults="./venturi_30/mesh" + wildcard(i) + "/" 

  call("mkdir -p " + dirResults, shell=True)
  dataFile = dirResults + "parameters_geo_october"
  call("cp parameters_geo_october " + dataFile, shell=True)

  setParameter("dirResults", dirResults, dataFile)
  setParameter("dirModel", "../../data/venturi/mor/mesh" + wildcard(i) + "/00061/", dataFile)

  if n_procs == 1:
    cmd = "./pbdw.exe " + dataFile
  else:
    cmd = "./pbdw.exe " + dataFile + " > " + dirResults + "or.log"

  print(cmd)
  call(cmd, shell=True)
  call("rm " + dirResults + "/*.geo", shell=True)
  call("rm " + dirResults + "/*.case", shell=True)
  call("rm " + dirResults + "/*.vct", shell=True)
  call("rm " + dirResults + "/*.scl", shell=True)

if n_procs == 1:
  for j in range(nbMeshes):
    loopSim(j)

else:
  p = Pool(processes = n_procs)
  p.map(loopSim, range(nbMeshes))
