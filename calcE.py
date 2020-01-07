import numpy as np
import numpy.linalg as LA
import numpy.matlib as matlib

def calcE_geometric(G, C, A, alpha):
    I = np.identity(A[:,0].size)
    return np.transpose(G) * alpha * LA.inv(I-(1-alpha)*A) * C * LA.inv(np.transpose(C)*C)

def calcE_geometric_noC(G, A, alpha):
    I = np.identity(A[:,0].size)
    return np.transpose(G) * alpha * LA.inv(I-(1-alpha)*A)


def calcE_geometric_noC_invA(G, invA, alpha):

    return np.transpose(G) * alpha * invA



def calcE_poisson():
    raise Exception("not implemented yet")
    return 

def get_invA(A, alpha):
    I = np.identity(A[:,0].size)
    return (I-(1-alpha)*A).I

def get_invApi(Api, alpha):
    Ipi = np.identity(Api[:,0].size)
    return LA.inv(Ipi-(1-alpha)*Api)

def get_invC(C):
    return LA.inv(np.transpose(C)*C)

def calcE_geometric_EP_invApiinvC(G, C, invC, invApi, alpha):
    return alpha * np.transpose(G) * C * invApi * invC



def calcE_geometric_EP(G, C, Api, alpha):
    Ipi = np.identity(Api[:,0].size)
#    print "G:", G[:,0], G[0,:]
#    print "C:", C[:,0], C[0,:]
#    print "C:", C[:,0], C[0,:]
    return alpha * np.transpose(G) * C * LA.inv(Ipi-(1-alpha)*Api) * LA.inv(np.transpose(C)*C)


def calcE_EPpoisson():
    raise Exception("not implemented yet")
    return


def main():
    C = np.matrix('1 0 ; 1 0 ; 1 0 ; 0 1 ; 0 1')
    Cplus = LA.inv(np.transpose(C)*C)
    G_0 = np.matrix('0 1 ; 1 0 ; 1 0 ; 0 1 ; 0 1')
    n = 5
    c = 2
    g = 2
    alpha = 0.1

    #Admits EP
    A_0 = np.matrix('.0 .4 .4 .1 .1 ; .4 .0 .4 .1 .1 ; .4 .4 .0 .1 .1 ; .1 .1 .1 .0 .7 ; .1 .1 .1 .7 .0')
    T_0 = np.matrix('1 .3 ; 0 .7')
    T_1 = np.matrix('.2 0 ; .8 1')

    I = np.identity(n)
    Ipi = np.identity(c)

    #Look at this. It's interpretable 
    Api = Cplus * np.transpose(C) * A_0 * C

    #should be equal and column stochastic
    print calcE_geometric(G_0, C, A_0, alpha)
    print calcE_EPgeometric(G_0, C, Api, alpha)
