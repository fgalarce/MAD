import os
import numpy as np
import matplotlib.pyplot as plt
from nibabel.testing import data_path

def writeVTKgrid(data, nx, ny, nz, name, aspect_ratio_x, aspect_ratio_y, aspect_ratio_z, origin_x, origin_y, origin_z):

  print("Writing " + name)
  f = open(name, "w")
  f.write("# vtk DataFile Version 2.0\n")
  f.write("Data animation\n")
  f.write("ASCII\n");
  f.write("DATASET STRUCTURED_POINTS\n");
  f.write("DIMENSIONS " + str(nx) + " " + str(ny) + " " + str(nz) + "\n");
  f.write("ASPECT_RATIO " + str(aspect_ratio_x) + " " + str(aspect_ratio_y) + " " + str(aspect_ratio_z) + "\n");
  f.write("ORIGIN " + str(origin_x) + " " + str(origin_y) + " " + str(origin_z) + "\n");
  f.write("POINT_DATA " + str(len(data)) + "\n");
  f.write("SCALARS Field float\n");
  f.write("LOOKUP_TABLE default\n");
  for h in data:
    if (abs(h) > 0.00000001):
      f.write(str(h) + "\n")
    else:
      f.write(str(0.0) + "\n")
  
  f.close()

#file = '/Users/galarce/research/mad/data/charite_data/Alfonso_20180424-142128/MPRAGE_0004/c1s2018-04-24_14-21-142923-00001-00176-1.nii'
file = '/Users/galarce/research/mad/data/charite_data/Alfonso_20180424-142128/MPRAGE_0004/c2s2018-04-24_14-21-142923-00001-00176-1.nii'

example_filename = os.path.join(data_path, file)

import nibabel as nib
img = nib.load(example_filename)
print("Raw data size")
print(img.shape)
img.get_data_dtype() == np.dtype(np.int16)

img.affine.shape
data = img.get_fdata()

##plt.imshow(data[:,:,100])
##plt.show()
x_min = -4.10216
x_max = 8.07514
y_min = -10.536
y_max = 6.3501
z_min = -6.68273
z_max = 6.72076

id_max_x = 0
id_max_y = 0
id_max_z = 0
id_min_x = 300
id_min_y = 300
id_min_z = 300
for k in range(np.shape(data)[2]):
  for j in range(np.shape(data)[1]):
    for i in range(np.shape(data)[0]):

      if (data[i,j,k] > 0.0 and i > id_max_x):
        id_max_x = i
      if (data[i,j,k] > 0.0 and j > id_max_y):
        id_max_y = j
      if (data[i,j,k] > 0.0 and k > id_max_z):
        id_max_z = k
      
      if (data[i,j,k] > 0.0 and i < id_min_x):
        id_min_x = i
      if (data[i,j,k] > 0.0 and j < id_min_y):
        id_min_y = j
      if (data[i,j,k] > 0.0 and k < id_min_z):
        id_min_z = k

print(id_min_x)
print(id_max_x)
print(id_min_y)
print(id_max_y)
print(id_min_z)
print(id_max_z)

data_short = np.zeros((id_max_y - id_min_y +1, id_max_x - id_min_x + 1, id_max_z - id_min_z + 1))

for k in range(id_min_z, id_max_z+1):
  for i in range(id_min_y, id_max_y+1):
    for j in range(id_min_x, id_max_x+1):
      data_short[i-id_min_y,j-id_min_x,k-id_min_z] = data[j, i, k]

#plt.imshow(data_short[:,:,60])
#plt.show()

scalex = (x_max - x_min) / float(np.shape(data_short)[0])
scaley = (y_max - y_min) / float(np.shape(data_short)[1])
scalez = (z_max - z_min) / float(np.shape(data_short)[2])

it = 0
print("Shrinked data size")
print(data_short.shape)
vtk_data = np.zeros(np.shape(data_short)[0]*np.shape(data_short)[1]*np.shape(data_short)[2] )
for k in range(np.shape(data_short)[2]):
  for j in range(np.shape(data_short)[1]):
    for i in range(np.shape(data_short)[0]):
      vtk_data[it] = data_short[i,j,k]
      it = it + 1

writeVTKgrid(vtk_data, np.shape(data_short)[0], np.shape(data_short)[1], np.shape(data_short)[2], "indicatrix.vtk", scalex, scaley, scalez, x_min, y_min, z_min)
