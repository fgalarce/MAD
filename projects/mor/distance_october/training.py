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

nbMeshes = 64

def loop(j):

  for i in range(nbMeshes):

    dirResults = "./training/"
    call("mkdir -p " + dirResults, shell=True)
    
    dataFile = dirResults  + "/parameters"
    
    call("cp parameters " + dataFile, shell=True)
    
    geometryData="../../data/venturi/mesh/venturi" + wildcard(i) + ".mesh"
    model="../../data/venturi/mor/mesh" + wildcard(j) + "/" + wildcard(i) + "/"
#    model="../../data/venturi/mor/mesh" + wildcard(j) + "/"
    maniFolder="../../data/venturi/dictionary/mesh" + wildcard(i) + "/"

    setParameter("model", model,  dataFile) 
    setParameter("geometryData", geometryData,  dataFile) 
    setParameter("maniFolder", maniFolder,  dataFile) 
    

    if n_procs == 1:
        cmd = "./training.exe " + dataFile + " " + str(j) + " " + str(i)
    else:
        cmd = "./training.exe " + dataFile + " " + str(j) + " " + str(i) + " > " + dirResults + "/distance_october.log"
   
    if (j != i): 
#    if (j == i): 
      print(cmd)
      call(cmd, shell=True) 
  
if n_procs == 1:
  for j in range(nbMeshes):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(nbMeshes)) 
