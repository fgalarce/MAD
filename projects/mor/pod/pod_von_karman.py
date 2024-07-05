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

n_procs = int(args.n_procs)

first=12
nbGeometries=20

mpi="/Home/flow/galarce/MAD/petsc/build/bin/mpirun"

def loop(i):

  maniFolder="../../../../data/von_karman/geo" + wildcard(i) + "/"
  geometryData="../../../../data/von_karman/fpc." + wildcard(i) + ".mesh"
  dirResults="../../../../data/von_karman/geo" + wildcard(i) + "/pod/"

  call("mkdir -p " + dirResults, shell=True)
  dataFile = dirResults + "parameters"
  call("cp par_street " + dataFile, shell=True)
  
  setParameter("maniFolder", maniFolder, dataFile) 
  setParameter("geometryData", geometryData, dataFile) 
  setParameter("dirResults", dirResults, dataFile) 
  
  if n_procs == 1:
      cmd = mpi + " -np 72 ./pod.exe " + dataFile
  else:
      cmd = mpi + " -np 12 ./pod.exe " + dataFile + " > " + dirResults + "/pod.log"
  
  print(cmd)
  call(cmd, shell=True) 

if n_procs == 1:
  for j in range(first,nbGeometries):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(first,nbGeometries)) 
