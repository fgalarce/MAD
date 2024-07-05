from subprocess import *
from numpy import *
from multiprocessing import Pool
import fileinput

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
nbSimulations=16

mpi="/Home/flow/galarce/MAD/source/external_software/petsc/build/bin/mpirun"

def loop(j):
  for n in range(2,105,5):
    dirResults="./pbdw_L2/test" + wildcard(j) + "/n" + str(n) + "/"
    dirField="../../pdes/elastodynamics/brain/target" + wildcard(j) + "/"
    call("mkdir -p " + dirResults, shell=True)
    dataFile = dirResults + "parameters"
    call("cp parameters_brain " + dataFile, shell=True)
    setParameter("dirResults", dirResults, dataFile) 
    setParameter("dirSyntheticField", dirField + "/disp" , dataFile) 
    setParameter("innerProduct", "H1" , dataFile) 
    setParameter("nbModes", n, dataFile) 
    nbIterations=getParameter(dirField + "/parameters", "nbIterations")
    setParameter("end", nbIterations, dataFile) 
    
    if n_procs == 1:
        cmd = mpi + " -np 72 ./pbdw.exe " + dataFile
    else:
        cmd = mpi + " -np 18 ./pbdw.exe " + dataFile + " > " + dirResults + "/elastodynamics.log"
    
    print(cmd)
    call(cmd, shell=True) 
 
if n_procs == 1:
  for j in range(nbSimulations):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(nbSimulations)) 
