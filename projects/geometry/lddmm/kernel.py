from pykeops.torch import Vi, Vj

def GaussKernel(sigma):
    x, y, b = Vi(0, 3), Vj(1, 3), Vj(2, 3)
    gamma = 1 / (sigma * sigma)
    D2 = x.sqdist(y)
    K = (-D2 * gamma).exp()
    return (K * b).sum_reduction(axis=1)

###################################################################
# Define "Gaussian-CauchyBinet" kernel :math:`(K(x,y,u,v)b)_i = \sum_j \exp(-\gamma\|x_i-y_j\|^2) \langle u_i,v_j\rangle^2 b_j`

def GaussLinKernel(sigma):
    x, y, u, v, b = Vi(0, 3), Vj(1, 3), Vi(2, 3), Vj(3, 3), Vj(4, 1)
    gamma = 1.0 / (sigma * sigma)
    D2 = x.sqdist(y)
    K = (-D2 * gamma).exp() * (u * v).sum() ** 2
    return (K * b).sum_reduction(axis=1)
