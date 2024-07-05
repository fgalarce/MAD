'''=============================================================================
  This file is part of the code MAD
  Multy-physics for biomedicAl engineering and Data assimilation.
  Copyright (C) 2017-2020,
    
     Felipe Galarce at INRIA

  MAD is free software; you can redistribute it and/or modify it under
  the terms of the GNU Lesser General Public License as published by the Free
  Software Foundation; either version 2.1 of the License, or (at your option)
  any later version.

  MAD is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
  more details.

  You should have received a copy of the GNU Lesser General Public License
  along with MAD. If not, see http://www.gnu.org/licenses/.
  ============================================================================='''

import numpy as np
import torch
import subprocess
import torch
import sys
from decimal import Decimal

def importEnsight(fileName, scaleFactor):

  print('.geo Ensight meshes are currently unavailable. Use .mesh instead')
  sys.exit()
  print('Importing Ensight file : ' + fileName + ".geo")

  # prepare loop
  fileIN = open(fileName + ".geo", "r")
  line = fileIN.readline()
  tria = np.array([0, 0, 0], dtype='int')
  countElemen = 0
  while line:
    if (line == 'coordinates\n'):
      nbPoints = int(fileIN.readline())
      vertices = np.zeros((nbPoints,3))
      for i in range(0,nbPoints):
        line = fileIN.readline()
        vertices[i,0] = float(line[0:11])
        vertices[i,1] = float(line[12:23])
        vertices[i,2] = float(line[24:35])

    if (line == 'tria3\n'): # This assumes a 3D mesh with surface made of triangles "tria3"
      nbElem = int(fileIN.readline())
      countElemen = countElemen + nbElem
      startIndex = tria.shape[0] - 1;
      tria.resize((nbElem + tria.shape[0], 3))
      for i in range(startIndex, nbElem + startIndex):
        line = fileIN.readline().split()
        tria[i,0] = int(line[0])
        tria[i,1] = int(line[1])
        tria[i,2] = int(line[2])

    line = fileIN.readline()

  # get only boundary elements
  pointsId = np.array([0])
  for i in range(0,tria.shape[0]):
    searchRepeated = np.where(pointsId == tria[i,0])
    if (searchRepeated[0].size == 0):
      pointsId = np.append(pointsId, tria[i,0])
    searchRepeated = np.where(pointsId == tria[i,1])
    if (searchRepeated[0].size == 0):
      pointsId = np.append(pointsId, tria[i,1])
    searchRepeated = np.where(pointsId == tria[i,2])
    if (searchRepeated[0].size == 0):
      pointsId = np.append(pointsId, tria[i,2])
  pointsId = np.delete(pointsId,0)

  # sort, set start from 0 index and export to Torch format
  pointsId = np.sort(pointsId) - 1
  v = vertices[pointsId,:]*scaleFactor;
  v = torch.from_numpy(v)
  t = torch.from_numpy(tria)
  torch.save(t, fileName + "_triangles.pt")
  torch.save(v, fileName + "_vertices.pt")
  subprocess.call("touch " + fileName + ".pt", shell=True)
  print("Written " + fileName+".pt")

def exportEnsight(fileName, q0, disp, triangles):

  print("Exporting displacement to Ensight.")

  # writing case 
  fileOut = open(fileName + ".case", 'w')
  fileOut.write("FORMAT\n")
  fileOut.write("type: ensight\n")
  fileOut.write("GEOMETRY\n")
  fileOut.write("model: 1 " + fileName + ".geo\n")
  fileOut.write("VARIABLE\n")
  fileOut.write("vector per node: 1 " + fileName + "2d " + fileName + "2d.*****.vct\n")
  fileOut.write("TIME\n")
  fileOut.write("time set: 1\n")
  fileOut.write("number of steps: 1\n")
  fileOut.write("filename start number: 0\n")
  fileOut.write("filename increment: 1\n")
  fileOut.write("time values:\n")
  fileOut.write("0\n")
  fileOut.close()

  # writing .vct
  nbVer = q0.shape[0]

  fileNameVct = fileName + ".00000.vct"

  fileOut = open(fileNameVct, "w")
  print("Writing "+fileNameVct)
  fileOut.write("Vector per node\n")

  for i in range(len(disp)):
    line = ""         
    if (disp[i][0] >= 0):
      line = line + " "
    line = line + str('%.5e' % Decimal(str(disp[i][0])))
    if ((i + 1) % 2 != 0):
      fileOut.write("\n")

  if (len(disp) % 2 != 0):
    fileOut.write(str('%.5e' % Decimal(str(disp[len(disp)-1][0]))) + " " + str('%.5e' % Decimal(str(disp[len(disp)-1][1]))) + " " + str('%.5e' % Decimal(str(disp[len(disp)-1][2]))) + "\n")

  fileOut.close()

  # writing .geo
  fileOut = open(fileName + ".geo", "w")
  fileOut.write("Geometry file\n")
  fileOut.write("Geometry file\n")
  fileOut.write("node id assign\n")
  fileOut.write("element id assign\n")
  fileOut.write("coordinates\n")
  fileOut.write(wildcard8(str(nbVer)) + "\n")
  for i in range(nbVer):
    line = ""
    for j in range(3):
      if (q0[i,j] >= 0):
        line = line + " "
      line = line + str('%.5e' % Decimal(str(q0[i,j])))
    line = line + "\n"
    fileOut.write(line)

  nbElem = triangles.shape[0]
  triangles = triangles + 1
  fileOut.write("part 1\nSurface\ntria3\n")
  fileOut.write(wildcard8(nbElem) + "\n")
  for i in range(nbElem):
    line = wildcard8(triangles[i,0]) + wildcard8(triangles[i,1]) + wildcard8(triangles[i,2]) + "\n"
    fileOut.write(line)
  fileOut.close()

def writeCase(folder, fileName, geoname, nt):

  fileNameCase = fileName + ".case"
  print("IO: Writing: " + fileNameCase)
  fileOut = open(folder + "/" +fileNameCase, 'w')
  fileOut.write("FORMAT\n")
  fileOut.write("type: ensight\n")
  fileOut.write("GEOMETRY\n")
  fileOut.write("model: 1 " + geoname + ".geo\n")
  fileOut.write("VARIABLE\n")
  fileOut.write("vector per node: 1 " + fileName + " " + fileName +".*****.vct\n")
  fileOut.write("TIME\n")
  fileOut.write("time set: 1\n")
  fileOut.write("number of steps: " + str(nt)  + "\n")
  fileOut.write("filename start number: 0\n")
  fileOut.write("filename increment: 1\n")
  fileOut.write("time values:\n")
  for i in range(nt):
    fileOut.write(str(i) + "\n")
  fileOut.close()

def writeGeo(folder, fileName, v_3d, t_3d, tetra, scaleFactor):

  v_3d = v_3d / scaleFactor
  nbVer = v_3d.shape[0]
  fileNameGeo = fileName + ".geo"
  print("IO: Writing : " + folder + fileNameGeo)
  fileOut = open(folder + "/" + fileName + ".geo", "w")
  fileOut.write("Geometry file\n")
  fileOut.write("Geometry file\n")
  fileOut.write("node id assign\n")
  fileOut.write("element id assign\n")
  fileOut.write("coordinates\n")
  fileOut.write("   " + str(nbVer) + "\n")
  for i in range(nbVer):
    line = str('%.5e' % Decimal(str(v_3d[i,0]))) + " " + str('%.5e' % Decimal(str(v_3d[i,1]))) + " " + str('%.5e' % Decimal(str(v_3d[i,2]))) + "\n"
    fileOut.write(line)

  nbElem = t_3d.shape[0]
  fileOut.write("\npart 1\nSurface\ntria3\n")
  fileOut.write("   " + str(nbElem) + "\n")
  for i in range(nbElem):
    line = str(t_3d[i,0]) + " " + str(t_3d[i,1]) + " " + str(t_3d[i,2]) + "\n"
    fileOut.write(line)

  nbElem = tetra.shape[0]
  fileOut.write("\npart 2\nVolume\ntetra4\n")
  fileOut.write("   " + str(nbElem) + "\n")
  for i in range(nbElem):
    line = str(tetra[i,0]) + " " + str(tetra[i,1]) + " " + str(tetra[i,2]) + " " + str(tetra[i,3]) + "\n"
    fileOut.write(line)

  fileOut.close()

def exportEnsight(folder, fileName, v_3d, t_3d, disp, v_map, tetra):

  print("IO: Exporting displacement to Ensight")

  # Writing case 
  fileNameCase = fileName + ".case"
  print("IO: Writing : " + folder + fileNameCase)
  fileOut = open(folder + "/" +fileNameCase, 'w')
  fileOut.write("FORMAT\n")
  fileOut.write("type: ensight\n")
  fileOut.write("GEOMETRY\n")
  fileOut.write("model: 1 " + fileName + ".geo\n")
  fileOut.write("VARIABLE\n")
  fileOut.write("vector per node: 1 " + fileName + " " + fileName +".*****.vct\n")
  fileOut.write("TIME\n")
  fileOut.write("time set: 1\n")
  fileOut.write("number of steps: 1\n")
  fileOut.write("filename start number: 0\n")
  fileOut.write("filename increment: 1\n")
  fileOut.write("time values:\n")
  fileOut.write("0\n")
  fileOut.close()

  # 2d to 3d
  disp3d = np.zeros(v_3d.shape)
  k = 0
  for i in range(v_map.shape[0]):
    disp3d[v_map[i]-1] = disp[i]
    
  # Writing .vct
  nbVer = v_3d.shape[0]
  fileNameVct = fileName + ".00000.vct"
  print("IO: Writing : " + folder + fileNameVct)
  fileOut = open(folder + "/" + fileNameVct, "w")
  fileOut.write("Vector per node\n")

  for i in range(len(disp3d)):
    line = ""
    for j in range(3):
      if (disp3d[i,j] >= 0):
        line = line + " "
      line = line + str('%.5e' % Decimal(disp3d[i,j]))
    
    if ((i + 1) % 2 != 0):
      line = line + ""
    else:
      line = line + "\n"

    fileOut.write(line)
  fileOut.close()

  # writing .geo
  fileNameGeo = fileName + ".geo"
  print("Writing : " + fileNameGeo)
  fileOut = open(folder + "/" + fileName + ".geo", "w")
  fileOut.write("Geometry file\n")
  fileOut.write("Geometry file\n")
  fileOut.write("node id assign\n")
  fileOut.write("element id assign\n")
  fileOut.write("coordinates\n")
  fileOut.write("   " + str(nbVer) + "\n")
  for i in range(nbVer):
    line = str('%.5e' % Decimal(str(v_3d[i,0]))) + " " + str('%.5e' % Decimal(str(v_3d[i,1]))) + " " + str('%.5e' % Decimal(str(v_3d[i,2]))) + "\n"
    fileOut.write(line)

  nbElem = t_3d.shape[0]
  fileOut.write("\npart 1\nSurface\ntria3\n")
  fileOut.write("   " + str(nbElem) + "\n")
  for i in range(nbElem):
    line = str(t_3d[i,0]) + " " + str(t_3d[i,1]) + " " + str(t_3d[i,2]) + "\n"
    fileOut.write(line)

  nbElem = tetra.shape[0]
  fileOut.write("\npart 2\nVolume\ntetra4\n")
  fileOut.write("   " + str(nbElem) + "\n")
  for i in range(nbElem):
    line = str(tetra[i,0]) + " " + str(tetra[i,1]) + " " + str(tetra[i,2]) + " " + str(tetra[i,3]) + "\n"
    fileOut.write(line)

  fileOut.close()

def writeState(folder, fileName, v_3d, disp, v_map, time):

  disp3d = np.zeros(v_3d.shape)
  k = 0
  for i in range(v_map.shape[0]):
    disp3d[v_map[i]-1] = disp[i]
    
  nbVer = v_3d.shape[0]
#  fileNameVct = fileName + "." + wildcard(time) + ".vct"
  fileNameVct = fileName + ".vct"
  print("IO: Writing: " + fileNameVct)
  fileOut = open(folder + "/" + fileNameVct, "w")
  fileOut.write("Vector per node\n")

  for i in range(len(disp3d)):
    line = ""
    for j in range(3):
      if (disp3d[i,j] >= 0):
        line = line + " "
      line = line + str('%.5e' % Decimal(disp3d[i,j]))
    
    if ((i + 1) % 2 != 0):
      line = line + ""
    else:
      line = line + "\n"

    fileOut.write(line)
  fileOut.close()
