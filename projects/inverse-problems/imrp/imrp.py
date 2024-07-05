from subprocess import *
from numpy import *
from multiprocessing import Pool
import fileinput

def fileReplace(top, bottom, fileName):
  tempFile =  open(fileName, 'r+' )
  for line in fileinput.input( fileName ):
      tempFile.write( line.replace(top, bottom)  )
  tempFile.close()

def getFromDataFile(dataFilePath, variableName):
  datafile = open(dataFilePath, 'r')
  line = datafile.readline()
  while line:
    if (line.split("=")[0] == variableName):
      return line.split("=")[1]
    line = datafile.readline()

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

def loop(j):

  dirResults='../../data/reconstruction/imrp/sim' + wildcard(j) + "/"
  call("mkdir -p " + dirResults, shell=True)
  
  dataFile = dirResults + "data"
  
  call("cp data " + dataFile, shell=True)
  
  fileReplace("dirResults=./results/", "dirResults=" + dirResults, dataFile) 
  
  if n_procs == 1:
      cmd = "./imrp.exe " + dataFile + " " + str(j)
  else:
      cmd = "./imrp.exe " + dataFile + " " + str(j) + " > " + dirResults + "/imrp.log"
  
  print("=======================================")
  print(cmd)
  print("=======================================")
  call(cmd, shell=True) 
  
if n_procs == 1:
  for j in range(16):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(16)) 
