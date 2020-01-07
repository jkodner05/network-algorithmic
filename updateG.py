import numpy as np
import numpy.linalg as LA
import numpy.random as nprandom
import numpy.matlib as matlib
#import numpy.matrix as matrix
from copy import copy
import random
import time

seed = "mustard"
random.seed(seed)
nprandom.seed(1234)

def Gassertions(C, commlangdistros, g):
    if commlangdistros[:,0].size != C[0].size:
        raise Exception("Wrong number of community language distributions")
    if (1-np.sum(commlangdistros)) > 0.000001:
        raise Exception("All speakers in all communities must speak something")

def initG_maxindicator(C, commlangdistros, innov=0):
    """learner in comm gets comm max likelihood lang
    with 1-innov probability
    else 
    jump to random lang
    with innov probability"""
    
    g = commlangdistros[0,:].size
    c = commlangdistros[:,0].size
    n = C[:,0].size
    Gassertions(C, commlangdistros, g)
    G = matlib.zeros((n,g))
    
    for comm in range(0,c):
        langdistros = commlangdistros[comm,:]
        maxg = np.argmax(langdistros)
        for h in range(0,C[:,0].size):
            if C[h,comm] == 1:
                if random.random() < innov:
                    G[h,random.randint(0,g-1)] = 1
                else:
                    G[h,maxg] = 1
    return G
    

def initG_indicator(hindices, commlangdistros, innov = 0):
    """individuals pick a single language
    distribute them across the population given the average distribution
    else
    individuals pick random language with innov probability"""

    g = commlangdistros[0,:].size
    c = commlangdistros[:,0].size
    n = hindices[0].shape[1]
#    print "n", hindices[0].shape[1]
#    Gassertions(C, commlangdistros, g)
    
    G = matlib.zeros((n,g))


    randmat = nprandom.rand(G.shape[0],1)

    for b in range(0,n):
        if randmat[b] < innov:
            G[b,random.randint(0,g-1)] = 1
        else:
            comm = np.asarray(hindices[1])[0][b]
            distr = np.asarray(commlangdistros[comm,:])[0]
            G[b,:] = nprandom.multinomial(1, distr)

    return G


def initG_indicator_old(C, commlangdistros, innov = 0):
    """individuals pick a single language
    distribute them across the population given the average distribution
    else
    individuals pick random language with innov probability"""

    g = commlangdistros[0,:].size
    c = commlangdistros[:,0].size
    n = C[:,0].size
    Gassertions(C, commlangdistros, g)
    
    G = matlib.zeros((n,g))

    for comm in range(0,c):
        lowerbound = 0
        rangedict = {}
        for i in range(0,(commlangdistros[comm,:].size)):
            prob = commlangdistros[comm,:].item(i)
            if prob == 0:
                continue
            rangedict[(lowerbound, lowerbound + prob)] = i
            lowerbound += prob
        for h in range(0,C[:,0].size):
            if C[h,comm] == 1:
                rando = random.random()
                for grange, grammar in rangedict.iteritems():
                    if grange[0] <= rando and rando < grange[1]:
                        if random.random() < innov:
                            G[h,random.randint(0,g-1)] = 1
                        else:
                            G[h,grammar] = 1
    return G

def initG_uniform(C, commlangdistros, innov = 0):
    """every individual internalizes grammars at rate equal to population rate
    with prob 1-innov
    else
    picks a random distribution NOT IMPLEMENTED"""

    g = commlangdistros[0,:].size
    c = commlangdistros[:,0].size
    n = C[:,0].size
    Gassertions(C, commlangdistros, g)
    
    G = matlib.zeros((n,g))
    eq1 = np.where(C == 1)

    G[eq1[0]] = commlangdistros[eq1[1]]

    return G


def initG_uniform_old(C, commlangdistros, innov = 0):
    """every individual internalizes grammars at rate equal to population rate
    with prob 1-innov
    else
    picks a random distribution NOT IMPLEMENTED"""

    g = commlangdistros[0,:].size
    c = commlangdistros[:,0].size
    n = C[:,0].size
    Gassertions(C, commlangdistros, g)
    
    G = matlib.zeros((n,g))
    for i in range(0,c):
        for h in range(0,C[:,0].size):
            if C[h,i] == 1:
                G[h] = commlangdistros[i,:]
    return G



def updateG_preservepop(C, G0, Gfunc, newcommlangdistros, prespopprob, innov=0):
    G1 = Gfunc(C, newcommlangdistros, innov)

    randmat = nprandom.rand(G1.shape[0])
    below = np.where(randmat < prespopprob)
    G1[below,:] = G0[below,:]

    return G1

def updateG_preservepop_old(C, G0, Gfunc, newcommlangdistros, prespopprob, innov=0):
    G1 = Gfunc(C, newcommlangdistros, innov)

    if C[0,:].size != C[:,0].size:
        for r in range(0,G0[:,0].size):
            if random.random() < prespopprob:
                G1[r] = G0[r]
    else:
        for r in range(0, C[:,0].size):
            if random.random() < prespopprob:
                G1[r] = G0[r]
    return G1

def updateG_movement(C, G0, Gfunc, newcommlangdistros, presprob, movements, innov=0):
    G1 = Gfunc(C, newcommlangdistros, innov)
    movedindices = []
    if C[0,:].size != C[:,0].size:
        for r in range(0,G0[:,0].size):
            moveprob = movements[r][0]
            moveddistro = movements[r][1]
            rand = random.random()
            if rand < moveprob:
                G1[r] = moveddistro
                movedindices.append(1)
            elif rand >= moveprob and rand < presprob+moveprob:
                G1[r] = G0[r]
                movedindices.append(0)
            else:
                movedindices.append(0)
    else:
        for r in range(0, C[:,0].size):
            moveprob = movements[r][0]
            moveddistro = movements[r][1]
            rand = random.random()
            if rand < moveprob:
                G1[r] = moveddistro
                movedindices.append(1)
            elif rand >= moveprob and rand < presprob+moveprob:
                G1[r] = G0[r]
                movedindices.append(0)
            else:
                movedindices.append(0)
    return G1, np.asarray(movedindices)


def main():
    n = 8
    commlangdistros = np.repeat(np.matrix([[0.0, 1.0]]), n, axis=0)
    C = createC(n, (1,)*n)
    G = initG_uniform(C, commlangdistros, 0) #initial lang distribution
    G1 = updateG_movement(C, G, initG_uniform, commlangdistros, 0, .5, np.matrix([[1.0,0.0]]))
    print G1
    exit()


    C = np.matrix('1 0 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 1 0 ; 0 1 ; 0 1')
    distros = np.matrix('1.0 0.0 ; 0.5 0.5')
#    C = np.matrix('1 0 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 1 0 ; 0 1 ; 0 1 ; 1 0 ; 1 0 ; 0 1 ; 0 1')
    print initG_indicator(C, distros)
    print ""
    print initG_uniform(C, distros)
    distros = np.matrix('1.0 0.0 ; 0.4 0.6')
    print initG_maxindicator(C, distros)


if __name__ == "__main__":
    main()
