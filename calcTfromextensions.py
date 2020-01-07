# tools for constructing transition matrices T

import numpy as np
import numpy.linalg as LA
import numpy.matlib as matlib

def paramhamming(params1, params2):
    return sum([1 for i, p1 in enumerate(params1) if p1 != params2[i]])

def createhammingdict(paramdict):
    hammingdict = {}
    for i, pi in paramdict.iteritems():
        for j, pj in paramdict.iteritems():
            hammingdict[(i,j)] = paramhamming(pi, pj)
    return hammingdict

def createT_neutral(numvars):
    return np.identity(numvars)

def createT_advantage_alpha(alpha):
    Ta = np.transpose(np.matrix([[1.0, 0.0], [1-alpha, alpha]]))
    return Ta
    
def createT_advantage_beta(beta):
    Tb = np.transpose(np.matrix([[beta, 1-beta], [0.0, 1.0]]))
    return Tb


def createT(target, paramdict, sentencedict, m=1):
    hammingdict = createhammingdict(paramdict)
    n = len(paramdict[1])
    normedsentencedict = {}
    for k, v in sentencedict.iteritems():
        normedsentencedict[k] = v.intersection(sentencedict[target])

    TtCalc = matlib.zeros((2**n, 2**n))

    for (i, j), dist in hammingdict.iteritems():
        TtCalc[i-1, j-1] = 0
        if dist == 1:
            if len(normedsentencedict[j].difference(normedsentencedict[i])) == 0:
                continue

            TtCalc[i-1,j-1] = 1.0/n*len(normedsentencedict[j].difference(normedsentencedict[i]))/len(sentencedict[target])

    for i, row in enumerate(TtCalc):
        TtCalc[i,i] = 1-np.sum(row)

    for i in range(0, m-1):
        TtCalc *= TtCalc

    #N & B use opposite convention, so we transpose
#    return TtCalc
    return np.transpose(TtCalc)
