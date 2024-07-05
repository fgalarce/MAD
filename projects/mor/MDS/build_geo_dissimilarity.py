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

import argparse
parser = argparse.ArgumentParser(description='parser')
parser.add_argument('-n', action="store", dest='n_procs', default=1)
args = parser.parse_args()

n_procs =int(args.n_procs)

nbMeshes = 4

def loop(j):

  for idMesh in range(nbMeshes):

    dirResults = "./results/"

    dataFile = dirResults + "/parameters" + str(idMesh) + "_" + str(j)
    
    call("mkdir -p " + dirResults, shell=True)
    call("cp parameters " + dataFile, shell=True)
    
    setParameter("surfaceMapping", "../../data/venturi/mappings/" + wildcard(j) + "/lddmm_" + wildcard(idMesh) + ".vct", dataFile) 
#    setParameter("surfaceMappingInverse", "../../data/venturi/mappings/" + wildcard(idMesh) + "/lddmm_" + wildcard(j) + ".vct", dataFile) 
    setParameter("surfaceMappingInverse", "../../data/venturi/mor/mesh" + wildcard(idMesh) + "/" + wildcard(j) + "/I_Tij." + wildcard(0) + ".vct", dataFile) 
    setParameter("geometryData", "../../data/venturi/mesh/venturi" + wildcard(j) + ".mesh",  dataFile) 
    setParameter("dirResults", dirResults,  dataFile) 

    if n_procs == 1:
        cmd = "./build_geo_dissimilarity.exe " + dataFile + " " + str(idMesh) + " " + str(j)
    else:
        cmd = "./build_geo_dissimilarity.exe " + dataFile + " " + str(idMesh) + " " + str(j) + " > " + "./L2_work_matrix/D" + str(idMesh) + "_" + str(j) + ".log"
   
    if (j != idMesh): 
      print(cmd)
      call(cmd, shell=True) 
  
if n_procs == 1:
  for j in range(nbMeshes):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(nbMeshes)) 
