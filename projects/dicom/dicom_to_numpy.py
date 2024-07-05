import vtk
from vtk.util import numpy_support
import numpy

# Define the parser
import argparse
parser = argparse.ArgumentParser(description='parser')
parser.add_argument('-f', action="store", dest='folderName', default=1)
parser.add_argument('-slice', action="store", dest='slice', default=1)
args = parser.parse_args()

PathDicom = args.folderName
reader = vtk.vtkDICOMImageReader()
reader.SetDirectoryName(PathDicom)
reader.Update()

# Load dimensions using `GetDataExtent`
_extent = reader.GetDataExtent()
ConstPixelDims = [_extent[1]-_extent[0]+1, _extent[3]-_extent[2]+1, _extent[5]-_extent[4]+1]

# Load spacing values
ConstPixelSpacing = reader.GetPixelSpacing()

# Get the 'vtkImageData' object from the reader
imageData = reader.GetOutput()
# Get the 'vtkPointData' object from the 'vtkImageData' object
pointData = imageData.GetPointData()
# Ensure that only one array exists within the 'vtkPointData' object
assert (pointData.GetNumberOfArrays()==1)
# Get the `vtkArray` (or whatever derived type) which is needed for the `numpy_support.vtk_to_numpy` function
arrayData = pointData.GetArray(0)

# Convert the `vtkArray` to a NumPy array
ArrayDicom = numpy_support.vtk_to_numpy(arrayData)
# Reshape the NumPy array to 3D using 'ConstPixelDims' as a 'shape'
ArrayDicom = ArrayDicom.reshape(ConstPixelDims, order='F')

print(ArrayDicom)
for i in range(ArrayDicom.shape[2]):
  numpy.savetxt("slice_" + str(i) + ".txt", ArrayDicom[0:int(ArrayDicom.shape[0]/3), int(2.0/3.0*ArrayDicom.shape[1]) : int(ArrayDicom.shape[1]), i]);
