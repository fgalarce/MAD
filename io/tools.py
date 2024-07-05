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

from subprocess import *
from numpy import *
from multiprocessing import Pool
import fileinput

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

def wildcard8(idIter):
  if (idIter < 10):
    iteration = "       " + str(idIter);
  elif (idIter < 100):
    iteration = "      " + str(idIter);
  elif (idIter < 1000):
    iteration = "     " + str(idIter);
  elif (idIter < 10000):
    iteration = "   " + str(idIter);
  elif (idIter < 100000):
    iteration = "  " + str(idIter);
  elif (idIter < 1000000):
    iteration = " " + str(idIter);
  elif (idIter < 10000000):
    iteration = str(idIter);
  return iteration
#
#def reverseNormals(lddmm, F, V, interior=1):
#  if dimension == 2:
#    V0, V1  = V.index_select(0, F[:, 0]), V.index_select(0, F[:, 1])
#    normals = zeros((size(V0,0), 3), dtype=float32);
#    for i in range(size(normals,0)):
#      normals[i][0] = - V0.detach().cpu().numpy()[i][1] + V1.detach().cpu().numpy()[i][1];
#      normals[i][1] = - V1.detach().cpu().numpy()[i][0] + V0.detach().cpu().numpy()[i][0];
#      normals[i][2] = 0.0;
#      normalization = sqrt(normals[i][0]*normals[i][0] + normals[i][1]*normals[i][1])
#      normals[i] = normals[i]/normalization;
#      N = 0.5 * (V0 + V1), torch.from_numpy(normals)
#  else:
#    V0, V1, V2 = V.index_select(0, F[:, 0]), V.index_select(0, F[:, 1]), V.index_select(0, F[:, 2])
#    N = .5 * torch.cross(V1 - V0, V2 - V0)
#  
#return normals;
#
