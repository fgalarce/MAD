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

nbSimulations=10
first=0
model_error=arange(0, 1e-4, 1e-4/nbSimulations)
kappa=random.uniform(1e-8, 1e-9, nbSimulations)
E=random.uniform(1e+5, 1e+6, nbSimulations)
nu=random.uniform(0.4, 0.45, nbSimulations)
#csf_factor=random.uniform(1.0, 1.1, nbSimulations)
#csf_factor=random.uniform(1.0, 1.02, nbSimulations)
csf_factor=random.uniform(1.08, 1.1, nbSimulations)
print("# simulations: " + str(nbSimulations))

def loop(j):

  dirResults="./model_error/sim" + wildcard(j) + "/"
  call("mkdir -p " + dirResults, shell=True)
  dataFile = dirResults + "parameters"
  call("cp par_brain " + dataFile, shell=True)

  setParameter("dirResults", dirResults, dataFile) 
  setParameter("inlet_u0", model_error[j], dataFile) 
#  setParameter("youngModulus", E[j], dataFile) 
#  setParameter("poissonRatio", nu[j], dataFile) 
#  setParameter("permeability", kappa[j], dataFile) 
#  setParameter("csf_factor", csf_factor[j], dataFile) 

  if n_procs == 1:
      cmd = mpi + " -np 72 ./poro.exe " + dataFile
  else:
      cmd = "./poro.exe " + dataFile + " > " + dirResults + "/poro.log"

  print(cmd)
  call(cmd, shell=True) 

if n_procs == 1:
  for j in range(first, nbSimulations):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(first, nbSimulations)) 
