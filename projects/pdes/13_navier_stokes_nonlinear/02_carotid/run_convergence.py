from subprocess import *
from numpy import *
from multiprocessing import Pool
import fileinput
import os
dir_results = os.environ['MAD_RESULTS']
dir_data = os.environ['MAD_DATA']

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

# Define the parser
import argparse
parser = argparse.ArgumentParser(description='parser')
parser.add_argument('-n', action="store", dest='n_procs', default=1)
args = parser.parse_args()

n_procs =int(args.n_procs)

timeStep=array([4e-3, 3e-3, 2e-3, 1e-3])
nbIterations=array([417, 556, 833, 1667])
meshSize=array([13,16,24,36,57,96])

mpirun="/home/server/home/01_research/MAD/petsc/build/bin/mpirun"
def loop(j):

  dirResults = dir_results + "/convergence/" + str(int(meshSize[j])) + "kv/"
  call("mkdir -p " + dirResults, shell=True)
  
  dataFile = dirResults + "data"
  
  call("cp par " + dataFile, shell=True)
  
  setParameter("dirResults", dirResults, dataFile) 
  setParameter("geometryData", dir_data + "/mesh/carotid/caroCritic_" + str(int(meshSize[j])) + "kv.mesh", dataFile) 
  
  if n_procs == 1:
      cmd = "./ns.exe " + dataFile
  else:
      cmd = mpirun + " -np 6 ./ns.exe " + dataFile + " > " + dirResults + "/ns.log"
  
  print(cmd)
  call(cmd, shell=True) 
  
if n_procs == 1:
  for j in range(len(meshSize)):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(len(meshSize))) 

