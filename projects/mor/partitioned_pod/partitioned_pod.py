from subprocess import *
from numpy import *
from multiprocessing import Pool
from fileinput import *

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
  for line in input( fileName ):
    if (line.split("=")[0] == parameter):
      datafile.write(line.replace(line[0:len(line)-1], line.split("=")[0] + "=" +  str(value)))
      succes = True
    else:
      datafile.write(line)
  datafile.close()
  if not succes:
    exit("ERROR: Parameter " + parameter + " not found in " + fileName)

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

mapping_n = loadtxt('../../data/model/pod_linear_H1/mapping_n.txt')

def loop(j):
  for i in range(5):

    dirResults='../../data/model/pod_linear_H1/window_time' + str(i) + "_HR" + str(j) + "/"
    call("mkdir -p " + dirResults, shell=True)
    
    dataFile = dirResults + "data"
    
    call("cp data " + dataFile, shell=True)
    
    setParameter("dirResults", dirResults, dataFile) 
#    setParameter("nbModes", mapping_n[j][i] + 50, dataFile) 
    
    if n_procs == 1:
        cmd = "./partitioned_pod.exe " + str(i) +" "+ str(j) + " " + dataFile
    else:
        cmd = "./partitioned_pod.exe " + str(i) +" "+ str(j) + " " + dataFile + " > " + dirResults + "/pod.log"
    
    print("=======================================")
    print(cmd)
    print("=======================================")
    call(cmd, shell=True) 
    
if n_procs == 1:
  for j in range(5):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(5)) 
