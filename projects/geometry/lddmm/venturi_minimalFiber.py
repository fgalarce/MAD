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
parser.add_argument('-np', action="store", dest='n_procs', default=1)
args = parser.parse_args()

n_procs =int(args.n_procs)
nbMeshes=64
nbTargets=8

def loop(idTarget):
  for j in range(nbMeshes):
    ensName="lddmm_target" + wildcard(j)
    dirResults="../../data/venturi/mappings/target" + wildcard(idTarget) + "/"
    call("mkdir -p " + dirResults, shell=True)
    source= "../../data/venturi/mesh/target" + wildcard(idTarget) + ".mesh"
    target="../../data/venturi/mesh/venturi" + wildcard(j) + ".mesh"

    run = "python3.7 LDDMM.py -scale=0.1 -save_folder=" + dirResults + " -name=" + ensName + " -source=" + source + " -target=" + target

    if n_procs == 1:
        cmd = run
    else:
        cmd = run + " > " + dirResults + "/lddmm_target" + wildcard(j) + ".log"
   
    print(cmd)
    call(cmd, shell=True)
  
if n_procs == 1:
  for j in range(nbTargets):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(nbTargets)) 
