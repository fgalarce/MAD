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
nbMeshes=71

def loop(idMeshTemp):
  for j in range(nbMeshes):

    dirResults="../../../../data/fpc/geo" + wildcard(j) + "/pod/piola/" + wildcard(idMeshTemp) + "/"
    templateGeometry="../../../../data/mesh/fpc/fpc." + wildcard(idMeshTemp) + ".mesh"
    geometryData="../../../../data/mesh/fpc/fpc." + wildcard(j) + ".mesh"
    dirSyntheticField="../../../../data/fpc/geo" + wildcard(idMeshTemp) + "/pod/velocity_mode"
    surfaceMapping="../../../../data/fpc/mappings/geo" + wildcard(idMeshTemp) + "/lddmm_" + wildcard(j) + ".00000.vct"
    surfaceMappingInverse="../../../../data/fpc/mappings/geo" + wildcard(j) + "/lddmm_" + wildcard(idMeshTemp) + ".00000.vct"

    dataFile = dirResults + "/par_interpolation"
    call("mkdir -p " + dirResults, shell=True)
    call("cp par_piola " + dataFile, shell=True)

    setParameter("dirResults", dirResults, dataFile) 
    setParameter("dirSyntheticField", dirSyntheticField, dataFile) 
    setParameter("geometryData", geometryData, dataFile) 
    setParameter("templateGeometry", templateGeometry, dataFile) 
    setParameter("surfaceMapping", surfaceMapping, dataFile) 
    setParameter("surfaceMappingInverse", surfaceMappingInverse, dataFile) 

    if n_procs == 1:
        cmd = mpi + " -np 72 ./piola.exe " + dataFile
    else:
        cmd = mpi + " -np 24 ./piola.exe " + dataFile + " > " + dirResults + "/piola.log"
   
    if (j != idMeshTemp): 
      print(cmd)
      call(cmd, shell=True) 
  
if n_procs == 1:
  for j in range(nbMeshes):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(nbMeshes)) 
