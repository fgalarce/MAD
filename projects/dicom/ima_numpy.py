import vtk
from vtk.util import numpy_support
import numpy

PathDicom = "/Users/galarce/research/mad/mad_data/charite_data/Vol9_IRMRI_IRMREsync/IRMREsync/"
PathDicom = "/Users/galarce/Downloads/Vol10/MPRAG/"
PathDicom = "/Users/galarce/research/mad/data/charite_data/Alfonso_20180424-142128/"
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
  numpy.savetxt("GW_slice_" + str(i) + ".txt", ArrayDicom[:,:,i]);

