from subprocess import *
from numpy import *
from multiprocessing import Pool
import fileinput
import os
dir_results = os.environ['MAD_RESULTS']

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

omega = 100
mesh = array([11,60,240])

def loop(j):

  dirResults=dir_results + "/biot_frecuency/delta/delta_1/layered_" +str(mesh[j]) + "kv/"
  call("mkdir -p " + dirResults, shell=True)
  
  dataFile = dirResults + "data"
  
  call("cp par " + dataFile, shell=True)
  
  setParameter("dirResults", dirResults, dataFile) 
  setParameter("geometryData", "MAD_DATA/mesh/layers/layer_" + str(mesh[j]) + "kv.mesh", dataFile) 
  setParameter("frecuency", str(omega), dataFile) 
  
  #mpirun="/home/server/home/02_research/MAD/petsc/build/bin/mpirun"
  mpirun="mpirun"
  if n_procs == 1:
      cmd = "mpirun -np 12 ./complex.exe " + dataFile
  else:
      cmd = mpirun + " -np 12 ./complex.exe " + dataFile + " > " + dirResults + "/complex.log"
  
  print("=======================================")
  print(cmd)
  print("=======================================")
  call(cmd, shell=True) 
  
if n_procs == 1:
  for j in range(len(mesh)):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(len(mesh))) 
