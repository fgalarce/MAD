from subprocess import *
import numpy as np
from multiprocessing import Pool
import fileinput
import sys

def importInriaMesh(fileNameInria):

  dimension = 3
  print("IO: Reading : " + fileNameInria)
  fileIN = open(fileNameInria, "r")
  line = fileIN.readline()
  while line:
    if len(line.split()) != 0:
      if line.split()[0] == 'Vertices':
        nbPoints = int(fileIN.readline())
        vertices = np.zeros((nbPoints,3))
        for i in range(0,nbPoints):
          line = fileIN.readline().split()
          vertices[i,0] = float(line[0])
          vertices[i,1] = float(line[1])
          vertices[i,2] = float(line[2]) 

        line = fileIN.readline()

    if len(line.split()) != 0:
      if line.split()[0] == 'Edges': 
        nbElem = int(fileIN.readline())
        edge = np.zeros((nbElem, 3), dtype=int)
        for i in range(0, nbElem):
          line = fileIN.readline().split()
          edge[i,0] = int(line[0])
          edge[i,1] = int(line[1])
          edge[i,2] = int(line[2])
    
        line = fileIN.readline()

    if len(line.split()) != 0:
      if line.split()[0] == 'Triangles': 
        nbElem = int(fileIN.readline())
        tria = np.zeros((nbElem, 4), dtype=int)
        for i in range(0, nbElem):
          line = fileIN.readline().split()
          tria[i,0] = int(line[0])
          tria[i,1] = int(line[1])
          tria[i,2] = int(line[2])
          tria[i,3] = int(line[3])
    
        line = fileIN.readline()

    if dimension == 3:

      if len(line.split()) != 0:
        if line.split()[0] == 'Tetrahedra': 
          nbElem = int(fileIN.readline())
          tetrahedra = np.zeros((nbElem, 5), dtype=int)
          for i in range(0, nbElem):
            line = fileIN.readline().split()
            tetrahedra[i,0] = int(line[0])
            tetrahedra[i,1] = int(line[1])
            tetrahedra[i,2] = int(line[2])
            tetrahedra[i,3] = int(line[3])
            tetrahedra[i,4] = int(line[4])
          line = fileIN.readline()
      
    line = fileIN.readline()

  fileIN.close()
  
  if dimension == 2:
    return vertices, edge, tria
#  return vertices, edge, tria, tetrahedra
  return vertices, tria, tetrahedra
