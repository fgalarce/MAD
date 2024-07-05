from subprocess import *
from numpy import *
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
      return line.split(" = ")[1]
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

noiseIter = 100
noise = arange(100,550,100)

# Only patients in HR3
nbSims = 4
sims = zeros((nbSims,), dtype = int)
sims[0] = 1
sims[1] = 3
sims[2] = 4
sims[3] = 13

def loopNoiseIter(j):
  for i in sims:
    for k in range(10,50,10):
      for noise_level in noise:

        """
        PBDW
        """
        targetFolder = "../../data/targets/sim" + wildcard(i) + "/"
        hr = float(getParameter(targetFolder + "data", "heartRate"))

        # Reconstruct first time window 
        UStimeStart=0.0
        UStimeEnd=0.2*60.0/hr
        dirResults="../../data/reconstruction/up_H1/noise/sim" + wildcard(i) + "/" + str(k) + "modes/cls_noise_" + str(int(noise_level)) + "_it" + str(j) + "/"
      
        call("mkdir -p " + dirResults, shell=True)
        dataFile = dirResults + "data"
        call("cp data " + dataFile, shell=True)
      
        setParameter("dirResults", dirResults, dataFile)
        setParameter("dirSyntheticField", targetFolder , dataFile)
        setParameter("dirMeasures2", "../../data/reconstruction/rr_H1/sim" + wildcard(i) + "/", dataFile)
        setParameter("UStimeStart", UStimeStart ,dataFile)
        setParameter("UStimeEnd", UStimeEnd ,dataFile)
        setParameter("heartRate", hr ,dataFile)
        setParameter("gaussianNoiseLevel", noise_level, dataFile)
        setParameter("nbModes", k, dataFile)
      
        if n_procs == 1:
          cmd = "./mixed_pbdw.exe " + dataFile
        else:
          cmd = "./mixed_pbdw.exe " + dataFile + " > " + dirResults + "or.log"

        print("it: " + str(j) + cmd)
        call(cmd, shell=True)

        """
        IMRP
        """
        call("cp data_pdrop " + dataFile, shell=True)
      
        setParameter("dirResults", dirResults, dataFile)
        setParameter("UStimeStart", UStimeStart ,dataFile)
        setParameter("UStimeEnd", UStimeEnd ,dataFile)

        if n_procs == 1:
          cmd = "./pdrop.exe " + dataFile
        else:
          cmd = "./pdrop.exe " + dataFile + " > " + dirResults + "pdrop.log"

        print("it: " + str(j) + cmd)
        call(cmd, shell=True)

        """
        clean up
        """
        call("rm " + dirResults + "/*.vct", shell=True)
        call("rm " + dirResults + "/*.scl", shell=True)
        call("rm " + dirResults + "/*.case", shell=True)
        call("rm " + dirResults + "/*.geo", shell=True)
    
if n_procs == 1:
  for j in range(noiseIter):
    loopNoiseIter(j)

else:
  p = Pool(processes = n_procs)
  p.map(loopNoiseIter, range(noiseIter))
