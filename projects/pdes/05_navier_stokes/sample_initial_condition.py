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

firstSim=0
firstMesh=0
nbSimulations=1
Re = arange(200,410,3.35)
r = arange(0.05, 0.1, 0.0007)
u_mean = 2.0 / 3.0 * 10;
nbMeshes=len(r)

print('# geometries: ' + str(nbMeshes))
print('# simulations / geometry : ' + str(len(Re)))

def loop(j):
  for idMesh in range(firstMesh, nbMeshes):

    geometryData="../../../../data/mesh/fpc/fpc." + wildcard(idMesh) + ".mesh"
    dirResults="../../../../data/fpc/geo" + wildcard(idMesh) + "/sim" + wildcard(j) 
    call("mkdir -p " + dirResults, shell=True)
    dataFile = dirResults + "/parameters"
    call("cp par_vonkarman " + dataFile, shell=True)
  
    setParameter("dirResults", dirResults, dataFile) 
    setParameter("geometryData", geometryData, dataFile) 
    setParameter("viscosity", 2 * r[idMesh] * u_mean / Re[j], dataFile) 

    if n_procs == 1:
        cmd = mpi + " -np 72 ./ns.exe " + dataFile
    else:
        cmd = "./ns.exe " + dataFile + " > " + dirResults + "/ns.log"
    
    print(cmd)
    call(cmd, shell=True) 

if n_procs == 1:
  for j in range(firstSim, nbSimulations):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(firstSim, nbSimulations)) 
