from subprocess import *
from numpy import *
from multiprocessing import Pool
import fileinput

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
first=0
nbSimulations=256

period = random.uniform(0.05, 0.1, nbSimulations)
poisson = random.uniform(0.4, 0.48, nbSimulations)
young = random.uniform(80000, 120000, nbSimulations)
amplitude = random.uniform(0.000002, 0.000005, nbSimulations)

mpi="/Users/galarce/research/mad/source/external_software/petsc/build/bin/mpirun"
mpi="/Home/flow/galarce/MAD/source/external_software/petsc/build/bin/mpirun"

def loop(j):

  dirResults="./brain/target" + wildcard(j) + "/"
  call("mkdir -p " + dirResults, shell=True)
  dataFile = dirResults + "parameters"
  call("cp parameters " + dataFile, shell=True)
  setParameter("dirResults", dirResults, dataFile) 

  setParameter("amplitude", amplitude[j], dataFile) 
  setParameter("poissonRatio", poisson[j], dataFile) 
  setParameter("youngModulus", young[j], dataFile) 
  setParameter("period", period[j], dataFile) 
  setParameter("nbIterations", int(period[j] / 0.001) + 1, dataFile) 
  
  if n_procs == 1:
      cmd = mpi + " -np 36 ./elasticity.exe " + dataFile
  else:
      cmd = "./elasticity.exe " + dataFile + " > " + dirResults + "/elastodynamics.log"
  
  print(cmd)
  call(cmd, shell=True) 

if n_procs == 1:
  for j in range(first, nbSimulations):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(nbSimulations)) 
