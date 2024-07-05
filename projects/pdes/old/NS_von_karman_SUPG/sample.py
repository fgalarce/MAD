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

def setParameterVector(parameter, value1, value2, fileName):
  datafile = open(fileName, 'r+')
  succes = False
  for line in fileinput.input( fileName ):
    if (line.split("=")[0] == parameter):
      datafile.write(line.replace(line[0:len(line)-1], line.split("=")[0] + "=" +  str(value1) + " " + str(value2) ))
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

mpi="/Home/flow/galarce/MAD/petsc/build/bin/mpirun"

Re = arange(50,202,2)
nbSimulations=len(Re)
r = arange(0.05, 0.092, 0.0006)
nbMeshes=70
meshesArray=array([34,33,32,31,30,29,28,24,23,22,21,20])

def loop(idMesh):
  for j in range(1,nbSimulations):

    geometryData="../../../../data/von_karman/fpc." + wildcard(idMesh) + ".mesh"
    dirResults="../../../../data/von_karman/geo" + wildcard(idMesh) + "/sim" + wildcard(j) 
    initial_condition_vel="../../../../data/von_karman/geo" + wildcard(idMesh) + "/sim" + wildcard(j - 1) + "/velocity.00149.vct"
    initial_condition_pre="../../../../data/von_karman/geo" + wildcard(idMesh) + "/sim" + wildcard(j - 1) + "/pressure.00149.scl"
    call("mkdir -p " + dirResults, shell=True)
    dataFile = dirResults + "/parameters"
    call("cp parameters " + dataFile, shell=True)
  
    setParameter("dirResults", dirResults, dataFile) 
    setParameter("initial_condition_vel", initial_condition_vel, dataFile) 
    setParameter("initial_condition_pre", initial_condition_pre, dataFile) 
    setParameter("geometryData", geometryData, dataFile) 
    setParameter("inlet_u0", Re[j]*0.001/(2*r[idMesh]), dataFile) 
    setParameter("nbIterations", 150, dataFile) 

    if n_procs == 1:
        cmd = mpi + " -np 72 ./ns.exe " + dataFile
    else:
        cmd = mpi + " -np 24 ./ns.exe " + dataFile + " > " + dirResults + "/ns.log"
    
    print(cmd)
    call(cmd, shell=True) 

if n_procs == 1:
  for j in meshesArray:
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, meshesArray) 
