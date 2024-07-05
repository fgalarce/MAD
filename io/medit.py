import numpy as np
import torch
import subprocess
import torch
import sys
from decimal import Decimal

def importMedit(fileName, scaleFactor, dimension, label=-1):
  
  # prepare loop
  fileNameMedit = fileName + ".mesh"
  print("IO: Reading : " + fileNameMedit)
  fileIN = open(fileNameMedit, "r")
  line = fileIN.readline()

  # edge_lddmm is the array that is going to be optimized
  edge_lddmm_id = np.array([0], dtype=int)

  while line:
    if len(line.split()) != 0:
      if line.split()[0] == 'Vertices':
        nbPoints = int(fileIN.readline())
        vertices = np.zeros((nbPoints,3))
        for i in range(0,nbPoints):
          line = fileIN.readline().split()
          vertices[i,0] = float(line[0]) * scaleFactor
          vertices[i,1] = float(line[1]) * scaleFactor
          vertices[i,2] = float(line[2]) * scaleFactor

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
          # filter elements to be optimized
          if (edge[i,2] == label):
            edge_lddmm_id = np.append(edge_lddmm_id, i)
    
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

  # get only boundary elements
  pointsId = np.array([0], dtype=int)
  for i in range(0,edge.shape[0]):
    searchRepeated = np.where(pointsId == edge[i,0])
    if (searchRepeated[0].size == 0):
      pointsId = np.append(pointsId, edge[i,0])
    searchRepeated = np.where(pointsId == edge[i,1])
    if (searchRepeated[0].size == 0):
      pointsId = np.append(pointsId, edge[i,1])

  pointsId_lddmm = np.array([0], dtype=int)
  edge_lddmm_id = np.delete(edge_lddmm_id,0)
  edge_lddmm = np.zeros((len(edge_lddmm_id),2), dtype=int)
  for i in range(len(edge_lddmm_id)):
    edge_lddmm[i][0] = edge[edge_lddmm_id[i]][0]
    edge_lddmm[i][1] = edge[edge_lddmm_id[i]][1]

  for i in range(0,edge_lddmm.shape[0]):
    searchRepeated = np.where(pointsId_lddmm == edge_lddmm[i,0])
    if (searchRepeated[0].size == 0):
      pointsId_lddmm = np.append(pointsId_lddmm, edge_lddmm[i,0])
    searchRepeated = np.where(pointsId_lddmm == edge_lddmm[i,1])
    if (searchRepeated[0].size == 0):
      pointsId_lddmm = np.append(pointsId_lddmm, edge_lddmm[i,1])

  pointsId_lddmm = np.delete(pointsId_lddmm,0)

  if dimension == 3:
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

  # sort, relabel and set start from 0 index 
  pointsId = np.sort(pointsId)
  pointsId_lddmm = np.sort(pointsId_lddmm)
  v = vertices[pointsId_lddmm-1,:]
  t = np.zeros(tria.shape)
  s = np.zeros(edge_lddmm.shape)

  for i in range(0,edge_lddmm.shape[0]):
    s[i,0] = (np.where(pointsId_lddmm == edge_lddmm[i,0]))[0]
    s[i,1] = (np.where(pointsId_lddmm == edge_lddmm[i,1]))[0]

  if dimension == 3:
    for i in range(0,tria.shape[0]):
      t[i,0] = (np.where(pointsId == tria[i,0]))[0]
      t[i,1] = (np.where(pointsId == tria[i,1]))[0]
      t[i,2] = (np.where(pointsId == tria[i,2]))[0]

  # keep 3d info
  v_3d = vertices
  t_3d = tria

  # export to torch (only those with the label)
  v = torch.from_numpy(v)
  s = torch.from_numpy(s)
  t = torch.from_numpy(t)

  if dimension == 2:
    return v, s, v_3d, edge, pointsId, pointsId_lddmm 
  else: 
    return v, t, v_3d, t_3d, pointsId, pointsId_lddmm
