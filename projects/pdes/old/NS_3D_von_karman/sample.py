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
    input("Parameter " + variableName + "not found, continue anyway?")
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
first=38

# average characteristic velocity
inletU0 = arange(0.5, 2.5, 0.02)
nbSimulations=len(inletU0)

print("Number of simulations : " + str(nbSimulations))

mpi="/Users/galarce/research/mad/source/external_software/petsc/build/bin/mpirun"
mpi="/Home/flow/galarce/MAD/source/external_software/petsc/build/bin/mpirun"

def loop(j):

  dirResults="/Home/flow/galarce/MAD/data/von_karman/geo00000/training_set/sim" + wildcard(j) + "/"
#
  call("mkdir -p " + dirResults, shell=True)
  dataFile = dirResults + "parameters"
  call("cp 2d_parameters " + dataFile, shell=True)
  setParameter("dirResults", dirResults, dataFile)

  setParameter("inlet_u0", inletU0[j], dataFile) 
  if (j > 0):
    setParameter("initial_condition_vel", "/Home/flow/galarce/MAD/data/von_karman/geo00000/training_set/sim" + wildcard(j-1) + "/velocity.00150.vct", dataFile)
    setParameter("initial_condition_pre", "/Home/flow/galarce/MAD/data/von_karman/geo00000/training_set/sim" + wildcard(j-1) + "/pressure.00150.scl", dataFile)
  
  if n_procs == 1:
      cmd = mpi + " -np 72 ./ns.exe " + dataFile
  else:
      cmd = "./ns.exe " + dataFile + " > " + dirResults + "/fpc.log"
  
  print(cmd)
  call(cmd, shell=True) 

if n_procs == 1:
  for j in range(first, nbSimulations):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(nbSimulations)) 
