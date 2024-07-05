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

n_procs =int(args.n_procs)

T = arange(0.0, 3, 0.25)
Bi = arange(0.0, 0.9, 0.1)

def loop(j):
  for B in Bi:

    dirResults="./resultsTclosure_" + str(int(100*T[j])) + "/Bi_" + str(int(B*10)) + "/"
    call("mkdir -p " + dirResults, shell=True)
    
    dataFile = dirResults + "data"
    
    call("cp par " + dataFile, shell=True)
    
    setParameter("dirResults", dirResults, dataFile) 
    setParameter("closureTime", T[j], dataFile) 
    setParameter("bingham", B, dataFile) 
    
    if n_procs == 1:
        cmd = "./hammer.exe " + dataFile
    else:
        cmd = "./hammer.exe " + dataFile + " > " + dirResults + "/log.log"
    
    print(cmd)
    call(cmd, shell=True) 
  
if n_procs == 1:
  for j in range(len(T)):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(len(T))) 
