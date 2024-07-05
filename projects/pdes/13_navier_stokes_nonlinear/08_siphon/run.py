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

nbsim=5
n = np.array(nbsim)
m = np.array(nbsim)
HR = np.array(nbSim)
HR[0]=35
n[0]=0.825350
m[0]=0.088821
HR[1]=45
n[1]=0.7754
m[1]=0.1482

def loop(j):

  dirResults="./results/HR" + str(HR[j]) + "/"
  call("mkdir -p " + dirResults, shell=True)
  
  dataFile = dirResults + "data"
  
  call("cp par " + dataFile, shell=True)
  
  setParameter("dirResults", dirResults, dataFile) 
  setParameter("power_law_m", m[j], dataFile) 
  setParameter("power_law_n", n[j], dataFile) 
  
  if n_procs == 1:
      cmd = "./ns.exe " + dataFile
  else:
      cmd = "./ns.exe " + dataFile + " > " + dirResults + "/<++>.log"
  
  print("=======================================")
  print(cmd)
  print("=======================================")
  call(cmd, shell=True) 
  
if n_procs == 1:
  for j in range(nbSim):
    loop(j)

else:
  p = Pool(processes = n_procs)
  p.map(loop, range(nbSim)) 
