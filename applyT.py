import numpy as np
import numpy.linalg as LA
import numpy.matlib as matlib

def calccombinedET(cE, Ts):
    T = np.zeros((Ts[0][0,:].size,Ts[0][0,:].size))
    for g in range(0,cE[:,0].size):
        T += cE[g,0] * Ts[g]

    return T

def cotcaught_p0():
    #cf Yang 2009
    cminus_a = 0.294
    cminus_b = 0.086
    cplus = 0.136
    p1 = 0.5
    p0 = ((-cplus*p1)/(p1-1) - cminus_b)/cminus_a
    return p0

def apply2varvariational_categorical_old(E, p0):
    G = matlib.zeros((E[0,:].size,E[:,0].size))

    for c in range(0,E[0,:].size):
        if E[0,c] > p0:
            G[c,0] = 1.0
        else:
            G[c,1] = 1.0

    return G


def apply2varvariational_categorical(E, p0):
    G = matlib.zeros((E.shape[1],E.shape[0]))

    above = np.where(E[0] > p0)
    notabove = np.where(E[0] <= p0)

    G[above[1],0] = 1
    G[notabove[1],1] = 1

    return G



def applyTneutral(E, Ts):
    return np.transpose(E)


def applyTinf(E, Ts):
    G = matlib.zeros((E[0,:].size,E[:,0].size))

    for c in range(0,E[0,:].size):
        T = calccombinedET(E[:,c],Ts)

        w, e = LA.eig(T)
        domeigv = np.where(abs(w-1) < 0.000001)[0][0]
#        print w, e

        Gc = e[:,domeigv]/np.sum(e[:,domeigv])
#        print c,":\t", np.imag(Gc)
        if np.sum(np.imag(Gc)) > 0.000000001:
            print "Got some imaginary stuff", Gc
            exit()
#        print G

        G[c,:] = Gc

    return G



def main():
    H = np.matrix('1 0 ; 1 0 ; 1 0 ; 0 1 ; 0 1')
    Hplus = LA.inv(np.transpose(H)*H)
    G_0 = np.matrix('0 1 ; 1 0 ; 1 0 ; 0 1 ; 0 1')
    n = 5
    c = 2
    g = 2
    alpha = 0.1

    A_0 = np.matrix('.0 .4 .4 .1 .1 ; .4 .0 .4 .1 .1 ; .4 .4 .0 .1 .1 ; .1 .1 .1 .0 .7 ; .1 .1 .1 .7 .0')
    T_0 = np.matrix('1 .3 ; 0 .7')
    T_1 = np.matrix('.2 0 ; .8 1')
    I = np.identity(n)
    Ipi = np.identity(c)

    #Gook at this. It's interpretable 
    Api = Hplus * np.transpose(H) * A_0 * H

    #should be equal and column stochastic
    E_1 = np.transpose(G_0) * alpha * LA.inv(I-(1-alpha)*A_0) * H * Hplus

    print applyTinf(E_1, (T_0, T_1))


