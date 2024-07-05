from subprocess import *
from sys import *
from multiprocessing import Pool
import fileinput
import argparse

parser = argparse.ArgumentParser(description='parser')
parser.add_argument('-n', action="store", dest='n_procs', default=1)
args = parser.parse_args()

n_procs = int(args.n_procs)
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

nbSims = 16

def loop(i):
  targetFolder = "../../data/targets/sim" + wildcard(i) + "/"

  hr = float(getParameter(targetFolder + "data", "heartRate"))

  UStimeStart=0.0
  UStimeEnd=60.0/hr

  dirResults = "../../data/reconstruction/u_H1/max_measures/sim" + wildcard(i) + "/"

  call("mkdir -p " + dirResults,shell=True)
  dataFile = dirResults + "data_noise"
  call("cp data " + dataFile, shell=True)

  setParameter("dirResults", dirResults , dataFile)
  setParameter("dirSyntheticField", targetFolder , dataFile)
  setParameter("UStimeStart", UStimeStart ,dataFile)
  setParameter("UStimeEnd", UStimeEnd ,dataFile)
  setParameter("heartRate", hr ,dataFile)

  if n_procs == 1:
    cmd = "./get_maxl.exe " + dataFile
  else:
    cmd = "./get_maxl.exe " + dataFile + " > " + dirResults + "/max_l.log"
  
  print("=======================================")
  print(cmd)
  print("=======================================")
  call(cmd, shell=True)
  call("rm " + dataFile, shell=True)

if n_procs == 1:
  for j in range(nbSims):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(nbSims)) 
