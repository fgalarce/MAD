from torch.autograd import grad

def RalstonIntegrator():
    def f(ODESystem, x0, nt, deltat=1.0):
        x = tuple(map(lambda x: x.clone(), x0))
        dt = deltat / nt
        l = [x]
        for i in range(nt):
            xdot = ODESystem(*x)
            xi = tuple(map(lambda x, xdot: x + (2 * dt / 3) * xdot, x, xdot))
            xdoti = ODESystem(*xi)
            x = tuple(map(lambda x, xdot, xdoti: x + (.25 * dt) * (xdot + 3 * xdoti), x, xdot, xdoti))
            l.append(x)
        return l
    return f

def Hamiltonian(K):
    def H(p, q):
        return .5 * (p * K(q, q, p)).sum()
    return H

def HamiltonianSystem(K):
    H = Hamiltonian(K)
    def HS(p, q):
        Gp, Gq = grad(H(p, q), (p, q), create_graph=True)
        return -Gq, Gp
    return HS
