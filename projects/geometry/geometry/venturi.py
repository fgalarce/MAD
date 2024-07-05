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

n_procs = int(args.n_procs)
nbMeshes=16

def loop(idMesh):
  for j in range(nbMeshes):
#  j = 1
    ensName="lddmm_" + wildcard(j)
    dirResults="../../data/venturi/mor/mesh" + wildcard(idMesh) + "/" + wildcard(j) + "/"
    call("mkdir -p " + dirResults, shell=True)

    # Template
    templateGeometry="../../data/venturi/mesh/venturi" + wildcard(idMesh) + ".mesh"
    surfaceMapping="../../data/venturi/mappings/" + wildcard(idMesh) + "/lddmm_" + wildcard(j) + ".vct"
    dirModel="../../data/venturi/mor/mesh" + wildcard(idMesh) + "/"

    # Target
    geometryData="../../data/venturi/mesh/venturi" + wildcard(j) + ".mesh"
    surfaceMappingInverse="../../data/venturi/mappings/" + wildcard(j) + "/lddmm_" + wildcard(idMesh) + ".vct"

    dataFile=dirResults + "/parameters" 
    call("cp parameters " + dataFile, shell=True)
    setParameter("dirResults", dirResults, dataFile) 
    setParameter("templateGeometry", templateGeometry, dataFile) 
    setParameter("surfaceMapping", surfaceMapping, dataFile) 
    setParameter("geometryData", geometryData, dataFile) 
    setParameter("surfaceMappingInverse", surfaceMappingInverse, dataFile) 
    setParameter("dirModel", dirModel, dataFile) 

    run = "./interpolate.exe " + dataFile

    if n_procs == 1:
        cmd = run
    else:
        cmd = run + " > " + dirResults + "/interp_piola." + wildcard(j) + ".log"
   
    if (j != idMesh): 
      print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
      print(cmd)
      print("- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ")
      call(cmd, shell=True) 
    
if n_procs == 1:
  for j in range(nbMeshes):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(nbMeshes)) 
