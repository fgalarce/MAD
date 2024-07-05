from numpy import *
from multiprocessing import Pool
from subprocess import *
import fileinput

def getFromDataFile(dataFilePath, variableName):
  datafile = open(dataFilePath, 'r')
  line = datafile.readline()
  while line:
    if (line.split("=")[0] == variableName):
      return line.split("=")[1]
    line = datafile.readline()

# search for line top and replace it by bottom
def fileReplace(top, bottom, fileName):
  tempFile =  open(fileName, 'r+' )
  for line in fileinput.input( fileName ):
      tempFile.write( line.replace(top, bottom)  )
  tempFile.close()

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

def test(i):

  timeStart = getFromDataFile("../../data/reconstruction/u_H1/sim" + wildcard(i) + "/data" , "UStimeStart")
  timeEnd = getFromDataFile("../../data/reconstruction/u_H1/sim" + wildcard(i) + "/data" , "UStimeEnd")
  dirResults = "../../data/vorticity/sim" + wildcard(i) + "/"
  call("mkdir -p " + dirResults, shell=True)
  call("cp data_carotid " + dirResults + "data", shell=True)

  fileReplace("dirResults=x", "dirResults=" + dirResults, dirResults + "data")
  fileReplace("UStimeStart=x", "UStimeStart=" + timeStart, dirResults + "data")
  fileReplace("UStimeEnd=x", "UStimeEnd=" + timeEnd, dirResults + "data")

  cmd = "./vorticity.exe " + dirResults + "data " + str(i) + " > " + dirResults + "/vorticity.log"
  print(cmd)
  call(cmd, shell=True)

p = Pool(processes = 4)
p.map(test, range(16))
