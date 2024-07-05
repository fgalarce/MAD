import matplotlib.pylab as plt
import scipy.sparse as sps
A = sps.rand(10000,10000, density=0.00001)
M = sps.csr_matrix(A)
plt.spy(M)
plt.show()

# Returns here '1.3.0'
matplotlib.__version__
