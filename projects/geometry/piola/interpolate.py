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
mpi="/Home/flow/galarce/MAD/petsc/build/bin/mpirun"

n_procs =int(args.n_procs)
nbMeshes=2

def loop(idMeshTemp):
  for j in range(nbMeshes):

    dirResults="../../../../data/von_karman/geo" + wildcard(idMeshTemp) + "/pod/" + wildcard(j) + "/"

    # [Template]
    templateGeometry="../../../../data/von_karman/fpc." + wildcard(idMeshTemp) + ".mesh"
    surfaceMapping="../../../../data/von_karman/mappings/" + wildcard(idMeshTemp) + "/lddmm_" + wildcard(j) + ".00000.vct"
    dirSyntheticField="../../../../data/von_karman/geo" + wildcard(idMeshTemp) + "/pod/mode"
    call("mkdir -p " + dirResults, shell=True)

    # [Target]
    geometryData="../../../../data/von_karman/fpc." + wildcard(j) + ".mesh"
    surfaceMappingInverse="../../../../data/von_karman/mappings/" + wildcard(j) + "/lddmm_" + wildcard(idMeshTemp) + ".00000.vct"

    dataFile = dirResults + "/parameters"
    call("cp par_interpolate " + dataFile, shell=True)

    setParameter("dirResults", dirResults, dataFile) 
    setParameter("templateGeometry", templateGeometry, dataFile) 
    setParameter("surfaceMapping", surfaceMapping, dataFile) 
    setParameter("dirSyntheticField", dirSyntheticField, dataFile) 
    setParameter("geometryData", geometryData, dataFile) 
    setParameter("surfaceMappingInverse", surfaceMappingInverse, dataFile) 

    if n_procs == 1:
        cmd = mpi + " -np 72 ./interpolate.exe " + dataFile
    else:
        cmd = mpi + " -np 24 ./interpolate.exe " + dataFile + " > " + dirResults + "/ns.log"
   
    if (j != idMeshTemp): 
      print(cmd)
      call(cmd, shell=True) 
  
if n_procs == 1:
  for j in range(nbMeshes):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(nbMeshes)) 
