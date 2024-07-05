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

n_procs=int(args.n_procs)

idMesh = 0
nbMeshes = 32

def loop(j):

  for n in range(1,20):
    dirResults="./results/n_" + str(n) + "/"
    call("mkdir -p " + dirResults, shell=True)
    
    dataFile = dirResults  + "/data"
    
    call("cp data " + dataFile, shell=True)
    
    templateModel="../../data/model/venturi/mesh" + wildcard(idMesh) + "/"
    geometryData="../../data/mesh/venturi/venturi" + wildcard(idMesh) + ".mesh"
    model="../../data/mesh/venturi/" + wildcard(j) + "/pod/" + wildcard(idMesh) + "/"
    mass_matrix_path="../../data/mesh/venturi/discretization_matrix/H1_mats" + wildcard(idMesh) + "/mass.bin"
    stiffness_matrix_path="../../data/mesh/venturi/discretization_matrix/H1_mats" + wildcard(idMesh) + "/stiff.bin"

    setParameter("dirResults", dirResults, dataFile) 
    setParameter("templateModel", templateModel,  dataFile) 
    setParameter("model", model,  dataFile) 
    setParameter("geometryData", geometryData,  dataFile) 
    setParameter("stiffness_matrix_path", stiffness_matrix_path,  dataFile) 
    setParameter("mass_matrix_path", mass_matrix_path,  dataFile) 
    setParameter("nbModes", n,  dataFile) 
    

    if n_procs == 1:
        cmd = "./grassmanian.exe " + dataFile + " " + str(j)
    else:
        cmd = "./grassmanian.exe " + dataFile + " " + str(j) + " > " + dirResults + "/haus" + str(idMesh) + "_" + str(n) + ".log"
   
    if (j != idMesh): 
      print("=======================================")
      print(cmd)
      print("=======================================")
      call(cmd, shell=True) 
  
if n_procs == 1:
  for j in range(1,nbMeshes):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(nbMeshes)) 
