import numpy as np
import numpy.linalg as LA
import numpy.matlib as matlib
import random
#from sinkhorn_knopp import sinkhorn_knopp as skp

seed = "mustard"
random.seed(seed)

def createA_onrouteoffroute_nonuniformcomms(C, n, csizes, sdfactor, kfactor, crosskfactor, numonroute, isthegreatdepression):
    """Models on-route / off-route distinction from Friedman 2012
    "Chicago" and "Midland" communities have outgoing connections only
    Chicago outgoing to on-route only
    Midlands outgoing to all
    On-route communities have outgoing to their four on-route neighbors and to one off-route
    Off-route have outgoing to their immediate off-route neighbors
    ----------
    randomly assigns connections 
    lower sdfactor = higher centrality
    higher kfactor = more intra-community connections
    higher crosskfactor = more inter-community connections
    number of communities on-route"""

    if sum(csizes) != n:
        raise Exception("Community sizes don't sum to n")

    A = matlib.zeros((n,n))

    #count connections
    baseindex = 0
    #community internal only
    #pick nodes to connect from normal distribution centered in center of clluster
    #sdfactor sets standard deviation
    for c, csize in enumerate(csizes):
        mu = csize/2
        sd = csize*sdfactor
        k = int(csize*kfactor)
        for edge in range(k):
            i = -1
            while i >= csize or i < 0:
                i = int(round(random.gauss(mu, sd)))
            j = i
            while j == i:
#                i = int(round(random.gauss(mu, sd)))
                j = random.randint(0,csize-1)
#            while j == i or j >= csize or j < 0:
#                j = int(round(random.gauss(mu, sd)))
            A[baseindex+i,baseindex+j] += 1
        baseindex += csize

    comms = np.nonzero(C)[1]
    
#    print A, "\n"
    
    #On route connections to Chicago
    if isthegreatdepression:
        baseindex = 0
        for c, csize in enumerate(csizes):
            if c == 0:
                baseindex += csize
                continue
            mu = csize/2
            sd = csize*sdfactor
            k = int(csize*crosskfactor)
            for edge in range(k):
                i = -1
                while i >= csize or i < 0:
                    i = random.randint(0,csize-1)
                j = i
                while j == i:
                    j = random.randint(0,csize-1)
                if c <= numonroute:
    #            if c == numonroute:
                    A[i,baseindex+j] += 1
                else:
                    A[sum(csizes[:-1])+i,baseindex+j] += 1
            baseindex += csize

    #on-route -> onroute connections
    baseindexi = 0
    for ci, cisize in enumerate(csizes):
        if ci > numonroute:
            break
        if ci == 0 or ci == len(csizes)-1:
            baseindexi += cisize
            continue
        baseindexj = 0
        for cj, cjsize in enumerate(csizes):
            if cj > numonroute:
                break
            if cj == 0 or cj == len(csizes)-1:
                baseindexj += cjsize
                continue
            if abs(cj-ci) < 2:
                baseindexj += cjsize
                continue
            mu = cisize/2
            sd = cisize*sdfactor
            k = int(csize*crosskfactor)
            for edge in range(k):
                i = -1
                while i >= cisize or i < 0:
                    i = random.randint(0,cisize-1)
                j = i
                while j == i:
                    j = random.randint(0,cjsize-1)
#                print baseindexi+i, baseindexj+j, "\t",  baseindexi, baseindexj
                A[baseindexi+i,baseindexj+j] += 1
            baseindexj += cjsize
#        print ""
        baseindexi += cisize

#    print A, "\n"

    #off-route -> off-route connections
    baseindexi = 0
    for ci, cisize in enumerate(csizes):
        if ci <= numonroute:
            baseindexi += cisize
            continue
        if ci == 0 or ci == len(csizes)-1:
            baseindexi += cisize
            continue
        baseindexj = 0
        for cj, cjsize in enumerate(csizes):
            if cj <= numonroute:
                baseindexj += cjsize
                continue
            if cj == 0 or cj == len(csizes)-1:
                baseindexj += cjsize
                continue
            if abs(cj-ci) < 2:
                baseindexj += cjsize
                continue
            mu = cisize/2
            sd = cisize*sdfactor
            k = int(csize*crosskfactor)
            for edge in range(k):
                i = -1
                while i >= cisize or i < 0:
                    i = random.randint(0,cisize-1)
                j = i
                while j == i:
                    j = random.randint(0,cjsize-1)
#                print baseindexi+i, baseindexj+j, "\t",  baseindexi, baseindexj
                A[baseindexi+i,baseindexj+j] += 1
            baseindexj += cjsize
#        print ""
        baseindexi += cisize

    #on-route -> off-route connections
    baseindexi = 0
    numoffroute = len(csizes)-2-numonroute
#    print "On->Off"
    for ci, cisize in enumerate(csizes):
        if ci > numonroute:
            break
        if ci == 0 or ci == len(csizes)-1:
            baseindexi += cisize
            continue
        connectindex = ci*numoffroute/numonroute+numonroute

        baseindexj = 0
        for cj, cjsize in enumerate(csizes):
            if cj <= numonroute:
                baseindexj += cjsize
                continue
            if cj == 0 or cj == len(csizes)-1:
                baseindexj += cjsize
                continue
            if cj != connectindex:
#            if abs(cj-connectindex) < 4:
                baseindexj += cjsize
                continue
#            print numonroute, ":", ci, cj
            mu = cisize/2
            sd = cisize*sdfactor
            k = int(csize*crosskfactor)*5
            for edge in range(k):
                i = -1
                while i >= cisize or i < 0:
                    i = random.randint(0,cisize-1)
                j = i
                while j == i:
                    j = random.randint(0,cjsize-1)
                #print baseindexi+i, baseindexj+j, "\t",  baseindexi, baseindexj
                if isthegreatdepression:
                    A[baseindexi+i,baseindexj+j] += 2
                else:
                    A[baseindexi+i,baseindexj+j] += 1
#                A[baseindexj+j,baseindexi+i] += 1
            baseindexj += cjsize
#        print ""
        baseindexi += cisize



#    print A
#    exit()

    A += np.identity(n)
#    print A
    #normalize columns to get percents
    colsums = np.sum(A,axis=0)

    A/= colsums
    
    return A

    

def createA_fullcorridor_nonuniformcomms(C, n, csizes, sdfactor, kfactor, crosskfactor, numonroute):
    """randomly assigns connections 
    lower sdfactor = higher centrality
    higher kfactor = more intra-community connections
    higher crosskfactor = more inter-community connections
    number of communities on-route"""

    if sum(csizes) != n:
        raise Exception("Community sizes don't sum to n")

    A = matlib.zeros((n,n))

    #count connections
    baseindex = 0
    #community internal only
    #pick nodes to connect from normal distribution centered in center of clluster
    #sdfactor sets standard deviation
    for c, csize in enumerate(csizes):
        mu = csize/2
        sd = csize*sdfactor
        k = int(csize*kfactor)
        for edge in range(k):
            i = -1
            while i >= csize or i < 0:
                i = int(round(random.gauss(mu, sd)))
            j = i
            while j == i:
#                i = int(round(random.gauss(mu, sd)))
                j = random.randint(0,csize-1)
#            while j == i or j >= csize or j < 0:
#                j = int(round(random.gauss(mu, sd)))
            A[baseindex+i,baseindex+j] += 1
        baseindex += csize

    comms = np.nonzero(C)[1]
    
#    print A, "\n"

    #On route connections to Chicago
    baseindex = 0
    for c, csize in enumerate(csizes):
        if c == 0:
            baseindex += csize
            continue
        mu = csize/2
        sd = csize*sdfactor
        k = int(csize*crosskfactor)
        for edge in range(k):
            i = -1
            while i >= csize or i < 0:
                i = random.randint(0,csize-1)
            j = i
            while j == i:
                j = random.randint(0,csize-1)
            if c <= numonroute:
                A[i,baseindex+j] += 1
            else:
                A[sum(csizes[:-1])+i,baseindex+j] += 1
        baseindex += csize

#    print A, "\n"

    #all other connections
    baseindexi = 0
    for ci, cisize in enumerate(csizes):
        if ci == 0 or ci == len(csizes)-1:
            baseindexi += cisize
            continue
        baseindexj = 0
        for cj, cjsize in enumerate(csizes):
            if cj == 0 or cj == len(csizes)-1:
                baseindexj += cjsize
                continue
            if abs(cj-ci) > 2:
                baseindexj += cjsize
                continue
            mu = cisize/2
            sd = cisize*sdfactor
            k = int(csize*crosskfactor)
            for edge in range(k):
                i = -1
                while i >= cisize or i < 0:
                    i = random.randint(0,cisize-1)
                j = i
                while j == i:
                    j = random.randint(0,cjsize-1)
#                print baseindexi+i, baseindexj+j, "\t",  baseindexi, baseindexj
                A[baseindexi+i,baseindexj+j] += 1
            baseindexj += cjsize
#        print ""
        baseindexi += cisize





    A += np.identity(n)
#    print A
    #normalize columns to get percents
    colsums = np.sum(A,axis=0)

    A/= colsums
    
    return A




def createA_merger_nonuniformcomms(C, n, csizes, sdfactor, kfactor, crosskfactor, numneighbors):
    """randomly assigns connections 
    lower sdfactor = higher centrality
    higher kfactor = more intra-community connections
    higher crosskfactor = more inter-community connections
    number of communities on-route"""

    if sum(csizes) != n:
        raise Exception("Community sizes don't sum to n")

    A = matlib.zeros((n,n))

    #count connections
    baseindex = 0
    #community internal only
    #pick nodes to connect from normal distribution centered in center of clluster
    #sdfactor sets standard deviation
    for c, csize in enumerate(csizes):
        mu = csize/2
        sd = csize*sdfactor
        k = int(csize*kfactor)
        for edge in range(k):
            i = -1
            while i >= csize or i < 0:
                i = int(round(random.gauss(mu, sd)))
            j = i
            while j == i:
#                i = int(round(random.gauss(mu, sd)))
                j = random.randint(0,csize-1)
#            while j == i or j >= csize or j < 0:
#                j = int(round(random.gauss(mu, sd)))
            A[baseindex+i,baseindex+j] += 1
        baseindex += csize

    comms = np.nonzero(C)[1]
    
#    print A, "\n"

    #On route connections to Innovation
    baseindex = 0
    for c, csize in enumerate(csizes):
        skip = random.random()
        if skip  < 1.0/2:
            print "skipping", c
            baseindex += csize
            continue
        mu = csize/2
        sd = csize*sdfactor
        k = int(csize*crosskfactor)
        for edge in range(k):
            i = -1
            while i >= csize or i < 0:
                i = random.randint(0,csize-1)
            j = i
            while j == i:
                j = random.randint(0,csize-1)
            A[sum(csizes[:-1])+i,baseindex+j] += 1
        baseindex += csize

#    print A, "\n"

    #all other connections
    baseindexi = 0
    for ci, cisize in enumerate(csizes):
        if ci == 0 or ci == len(csizes)-1:
            baseindexi += cisize
            continue
        baseindexj = 0
        neighbors = random.sample(range(0,len(csizes)),numneighbors)
        for cj, cjsize in enumerate(csizes):
            if cj == 0 or cj == len(csizes)-1:
                baseindexj += cjsize
                continue
            if cj not in neighbors:
                baseindexj += cjsize
                continue
#            if abs(cj-ci) > numneighbors:
#                baseindexj += cjsize
#                continue
            mu = cisize/2
            sd = cisize*sdfactor
            k = int(csize*crosskfactor)
#            if random.random() > 1/numneighbors:
#                k /= numneighbors
#                baseindexj += cjsize
#                continue
            for edge in range(k):
                i = -1
                while i >= cisize or i < 0:
                    i = random.randint(0,cisize-1)
                j = i
                while j == i:
                    j = random.randint(0,cjsize-1)
#                print baseindexi+i, baseindexj+j, "\t",  baseindexi, baseindexj
                A[baseindexi+i,baseindexj+j] += 1
            baseindexj += cjsize
#        print ""
        baseindexi += cisize





    A += np.identity(n)
#    print A
    #normalize columns to get percents
    colsums = np.sum(A,axis=0)

    A/= colsums
    
    return A






def createA_simplecorridor_nonuniformcomms(C, n, csizes, sdfactor, kfactor, stage1CCIcrossk, stage2CCIcrossk, MIDcrossk):
    """randomly assigns connections 
    lower sdfactor = higher centrality
    higher kfactor = more intra-community connections
    higher crosskfactor = more inter-community connections
    higher bifactor = more bidirectionality, 1 = non-directed graph"""

    if sum(csizes) != n:
        raise Exception("Community sizes don't sum to n")

    A1 = matlib.zeros((n,n))
    A2 = matlib.zeros((n,n))

    #count connections
    baseindex = 0
    #community internal only
    #pick nodes to connect from normal distribution centered in center of clluster
    #sdfactor sets standard deviation
    for c, csize in enumerate(csizes):
        mu = csize/2
        sd = csize*sdfactor
        k = int(csize*kfactor)
        print "cluster", c, "internal:\t", k
        for edge in range(k):
            i = -1
            while i >= csize or i < 0:
                i = int(round(random.gauss(mu, sd)))
            j = i
            while j == i:
                j = random.randint(0,csize-1)
#            while j == i or j >= csize or j < 0:
#                j = int(round(random.gauss(mu, sd)))
            A1[baseindex+i,baseindex+j] += 1
            A2[baseindex+i,baseindex+j] += 1
        baseindex += csize

    comms = np.nonzero(C)[1]
    
    A1 += np.identity(n)
    A2 += np.identity(n)

    #community external
    #MID
    k = int(MIDcrossk*csizes[1])
    print "MID -> Corridor:\t", k
    for edge in range(k):
        i = -1
        while i >= n or i < 0:
            i = random.randint(csizes[0]+csizes[1],n-1)
        j = -1
        while comms[:,i] == comms[:,j] or j >=  csizes[0]+csizes[1] or j < csizes[0]:
            j = random.randint(csizes[0],csizes[0]+csizes[1]-1)
        A1[i,j] += int(csizes[1]*kfactor)
        A2[i,j] += int(csizes[1]*kfactor)



    k = int(stage1CCIcrossk*csizes[1])
    print "CCI 1 -> Corridor:\t", k
    for edge in range(k):
        i = -1
        while i >= n or i < 0:
            i = random.randint(0,csizes[0]-1)
        j = -1
        while comms[:,i] == comms[:,j] or j >=  csizes[0]+csizes[1] or j < csizes[0]:
            j = random.randint(csizes[0],csizes[0]+csizes[1]-1)
        A1[i,j] += int(csizes[1]*kfactor)


    k = int(stage2CCIcrossk*csizes[1])
    print "CCI 2 -> Corridor:\t", k
    for edge in range(k):
        i = -1
        while i >= n or i < 0:
            i = random.randint(0,csizes[0]-1)
        j = -1
        while comms[:,i] == comms[:,j] or j >=  csizes[0]+csizes[1] or j < csizes[0]:
            j = random.randint(csizes[0],csizes[0]+csizes[1]-1)
        A2[i,j] += int(csizes[1]*kfactor)

    
#    A1 = np.transpose(A1)
#    A2 = np.transpose(A2)
    print A1
    print A2
    #normalize columns to get percents
    colsums = np.sum(A1,axis=0)
    A1/= colsums
    colsums = np.sum(A2,axis=0)
    A2/= colsums
    
    return A1, A2



def createA_simplecorridor_uniform(n, csizes, stage1CCI, stage2CCI, MID):
    """comm 1 chicago, comm 2 corridor, comm 3 midlands"""

    if sum(csizes) != n:
        raise Exception("Community sizes don't sum to n")
    elif len(csizes) != 3:
        raise Exception("There need to be three communities")
        
    A1 = matlib.zeros((n,n))
    A2 = matlib.zeros((n,n))

    #Corridor
    internal1 = (1.0-stage1CCI-MID)/(csizes[1])
    internal2 = (1.0-stage2CCI-MID)/(csizes[1])
    fromCCI1 = stage1CCI/csizes[1]
    fromCCI2 = stage2CCI/csizes[1]
    fromMID = MID/csizes[1]
    for i in range(0,n):
        for j in range(csizes[1],csizes[1]+csizes[2]):
            if i >= csizes[1] and i < csizes[1]+csizes[2]:
                A1[i,j] = internal1
                A2[i,j] = internal2
            elif i < csizes[1]: #Chicago
                A1[i,j] = fromCCI1
                A2[i,j] = fromCCI2
            else: #Midlands
                A1[i,j] = fromMID
                A2[i,j] = fromMID

    #Chicago
    internalCCI = 1.0/csizes[0]
    for i in range(0,csizes[0]):
        for j in range(0,csizes[0]):
            A1[i,j] = internalCCI
            A2[i,j] = internalCCI

    #Midlands
    internalMID = 1.0/csizes[0]
    for i in range(n-csizes[2],n):
        for j in range(n-csizes[2],n):
            A1[i,j] = internalMID
            A2[i,j] = internalMID

    return A1, A2


def createA_weaklyconnectedcomms(n, csizes, numouts, intraweight, sdfactor, kfactor):
    """like createA_nonuniformcommunities,
    but attempts k*csize connections per community 
    uniform weights. Will abort early if sdfactor is too small to make that happen
    attempt normalization with Sinkhorn-Knopp algorithm
    undirected if bifactor = 1"""

    if sum(csizes) != n:
        raise Exception("Community sizes don't sum to n")

    A = np.identity(n)

    internal = intraweight
    external = 1.0-intraweight

    #count connections
    baseindex = 0
    #community internal only
    #pick nodes to connect from normal distribution centered in center of clluster
    #sdfactor sets standard deviation
    for c, csize in enumerate(csizes):
        mu = csize/2
        sd = csize*sdfactor
        k = int(csize*kfactor)
        for edge in range(k):
            i = -1
            j = -1
            tries = 0
#            print edge
            while A[baseindex+i,baseindex+j] > 0 or (i < 0 and j < 0):
                #abort once a community is fully connected
                tries += 1
                if edge+1 > (csize**2-csize)/2 or tries > csize**2:
                    toomany = True
                    break
                i = -1
                while i >= csize or i < 0:
                    i = int(round(random.gauss(mu, sd)))
                j = -1
                while j == -1:
                    j = random.randint(0,csize-1)
            A[baseindex+i,baseindex+j] = internal

            i = -1
            j = -1
            while A[baseindex+i,baseindex+j] > 0 or (i < 0 and j < 0):
                #abort once a community is fully connected
                tries += 1
                if edge+1 > (csize**2-csize)/2 or tries > csize**2:
                    toomany = True
                    break
                j = -1
                while j >= csize or j < 0:
                    j = int(round(random.gauss(mu, sd)))
                i = -1
                while i == -1:
                    i = random.randint(0,csize-1)
            A[baseindex+i,baseindex+j] = internal

        for edge in range(0, numouts):
            i = random.randint(baseindex,baseindex+csize-1)
            j = baseindex
            while j >= baseindex and j < baseindex+csize:
                j = random.randint(csizes[0],n-1)
            A[i,j] = external

        baseindex += csize


    print A
    print A[0:csizes[0],0:csizes[0]]
    colsums = np.sum(A,axis=0)
    A/= colsums

    return A


def createA_adopters_continuous(n, adoptercatgrowth, extrafactor = 0):
    """Creates a network with adopter categories (early adopters...laggards) cf Rogers 1995
    speakers consider input only from those near and above them
    they speak to those immediately below them
    speakers with low indices are innovators and speakers with high indices are laggards
    adoptercatgrowth defines how big the categories grow as things go along
    0 no growth
    extrafactor*n how many extra random connections there are

    The difference here is that there is a continuous categorization rather than a categorical one"""

#    A = matlib.zeros((n,n))
    A = np.identity(n)
    basesd = 0.01 
    for i in range(0, n-1):
        mu = i+1
        sd = basesd+basesd*(adoptercatgrowth*i*n)
        for k in range(0,n):
            j = i
            while j < 0 or j == i or j > n-1:
                j = int(round(random.gauss(mu, sd)))
                if j <= i:
                    j = i + (i-j)
            A[i,j] += 1   
            A[i,i] == 1

    for k in range(0,int(extrafactor*n)):
        i = 0
        j = 0
        while i == 0:
            i = random.randint(0,n-1)        
        while j == 0:
            j = random.randint(0,n-1)        
        A[i,j] += 1
        
    print A
    colsums = np.sum(A,axis=0)
    A/= colsums

    return A

def createA_nonuniformK(C, n, csizes, sdfactor, kfactor, crosskfactor, bifactor=1.0):
    """like createA_nonuniformcommunities,
    but attempts k*csize connections per community 
    uniform weights. Will abort early if sdfactor is too small to make that happen
    attempt normalization with Sinkhorn-Knopp algorithm
    undirected if bifactor = 1"""

    if sum(csizes) != n:
        raise Exception("Community sizes don't sum to n")

    A = matlib.zeros((n,n))

    #count connections
    baseindex = 0
    #community internal only
    #pick nodes to connect from normal distribution centered in center of clluster
    #sdfactor sets standard deviation
    for c, csize in enumerate(csizes):
        mu = csize/2
        sd = csize*sdfactor
        k = int(csize*kfactor)
        for edge in range(k):
            i = 0
            j = 0
            tries = 0
            toomany = False
#            print edge
            while A[baseindex+i,baseindex+j] == 1 or A[baseindex+j,baseindex+i] == 1 or i == j:
                #abort once a community is fully connected
                tries += 1
                if edge+1 > (csize**2-csize)/2 or tries > csize**2:
                    toomany = True
                    break
                i = -1
                while i >= csize or i < 0:
                    i = int(round(random.gauss(mu, sd)))
                j = i
                while j == i:
                    j = random.randint(0,csize-1)
            if not toomany:
                A[baseindex+i,baseindex+j] += 1
        baseindex += csize

    comms = np.nonzero(C)[1]
#    print communitys
    #community external
    #pick nodes uniformly
    #only consider pairs that are in separate communitys
    openslots = (n**2-sum([csize**2 for csize in csizes]))/2
    if len(csizes) > 1:
        k = int(crosskfactor*n)
        for edge in range(k):
            if edge+1 > openslots:
                break
            i = 0
            j = 0
            while A[i,j] == 1 or A[j,i] == 1 or i == j:
                i = -1
                while i >= n or i < 0:
                    i = random.randint(0,n-1)
                j = -1
                while comms[:,i] == comms[:,j] or j >=  n or j < 0:
                    j = random.randint(0,n-1)
            A[i,j] += 1

    #undirected graph
    A += (np.transpose(A)*bifactor)


    #normalize
#    sk = skp.SinkhornKnopp()
#    A = sk.fit(A)

    for i in range(0, 5*n):
        rowsums = np.sum(A,axis=1)
        A /= rowsums
        colsums = np.sum(A,axis=0)
        A/= colsums
    return A

    
def createA_nonuniformcomms(C, n, csizes, sdfactor, kfactor, crosskfactor, bifactor):
    """randomly assigns connections 
    lower sdfactor = higher centrality
    higher kfactor = more intra-community connections
    higher crosskfactor = more inter-community connections
    higher bifactor = more bidirectionality, 1 = non-directed graph"""

    print "USING UNIFORMLY DISTRIBUTED INTRACLUSTER CONNECTIONS!"

    if sum(csizes) != n:
        raise Exception("Community sizes don't sum to n")

    A = matlib.zeros((n,n))

    #count connections
    baseindex = 0
    #community internal only
    #pick nodes to connect from normal distribution centered in center of clluster
    #sdfactor sets standard deviation
    for c, csize in enumerate(csizes):
        mu = csize/2
        sd = csize*sdfactor
        k = int(csize*kfactor)
        for edge in range(k):
            i = -1
            while i >= csize or i < 0:
                i = random.randint(0,csize-1)
#                i = int(round(random.gauss(mu, sd)))
            j = i
            while j == i:
                j = random.randint(0,csize-1)
#            while j == i or j >= csize or j < 0:
#                j = int(round(random.gauss(mu, sd)))
            A[baseindex+i,baseindex+j] += 1
        baseindex += csize

    comms = np.nonzero(C)[1]


    #community external
    #pick nodes uniformly
    #only consider pairs that are in separate communities
    if len(csizes) > 1:
        k = int(crosskfactor*n)
        for edge in range(k):
            i = -1
            while i >= n or i < 0:
                i = random.randint(0,n-1)
            j = -1
            while comms[i] == comms[j] or j >=  n or j < 0:
                j = random.randint(0,n-1)
            A[i,j] += 1

    #set diagonals to 1 so we can normalize safely
#    print np.sum(A,axis=1)

#    A = np.transpose(A)
    A += bifactor*np.transpose(A)
    A += np.identity(n)
    print A
    #normalize columns to get percents
    colsums = np.sum(A,axis=0)

    A/= colsums
    
    return A


def createA_uniformcomms(n, csizes, crosscprobs):
    """uniform within community, uniform between community distributed across crosscprobs probability masses"""
    if sum(csizes) != n:
        raise Exception("Community sizes don't sum to n")

#    if crosscprob < 0 or crosscprob > 1:
#        raise Exception("Inter-community jump prob needs to be a probability")
       
    A = matlib.zeros((n,n))

    baseindex = 0
    for c, csize in enumerate(csizes):
        if len(csizes) > 1:
            internal = (1.0-crosscprobs[c])/(csize)
            external = (crosscprobs[c])/(n-csize)
        else:
            internal = 1.0/csize
            external = 0
        for i in range(0,n):
            for j in range(baseindex,csize+baseindex):
                if i >= baseindex and i < csize+baseindex:
                    A[i,j] = internal
#                    if i != j:
#                        A[i,j] = internal
                else:
                    A[i,j] = external
        baseindex += csize

    return A

def createA_weakandstrong(n, csizes, crosscprobs,weakpercent):
    """like createA_uniformcommunitys except only weakpercent of the nodes in a community have cross-community connections"""
    if sum(csizes) != n:
        raise Exception("Community sizes don't sum to n")

    A = matlib.zeros((n,n))

    baseindex = 0
    for c, csize in enumerate(csizes):
        internalstrong = 1.0/csize
        externalstrong = 0
        internalweak = internalstrong
        externalweak = externalstrong
        if len(csizes) > 1:
            internalweak = (1.0-crosscprobs[c])/(csize)
            externalweak = crosscprobs[c]/(n-csize)
        for i in range(0,n):
            for j in range(baseindex,csize+baseindex):
                if i >= baseindex and i < csize+baseindex:
                    if j-baseindex < weakpercent*csize:
                        A[i,j] = internalweak
                    else:
                        A[i,j] = internalstrong
#                    if i != j:
#                        A[i,j] = internal
                else:
                    if j-baseindex < weakpercent*csize:
                        A[i,j] = externalweak
                    else:
                        A[i,j] = externalstrong

        baseindex += csize

    return A


def createC(n, csizes):
    if 0 in csizes:
        raise Exception("Can't have 0-size communitys")
    if len(csizes) == 0:
        raise Exception("Must have at least one community")
    if sum(csizes) != n:
        raise Exception("Community sizes don't sum to n")
        
    c = len(csizes)
    C = matlib.zeros((n,c))
    baseindex = 0
    for j, csize in enumerate(csizes):
        C[baseindex:baseindex+csize,j] = 1.0
        baseindex += csize
    return C

def createApi(A,C):
    return LA.inv(np.transpose(C)*C) * np.transpose(C) * A * C


def main():

    csizes = (4,4,4,4)
    n = sum(csizes)
    C = createC(n, csizes)
    print "A1\n"
    A1 = createA_onrouteoffroute_nonuniformcomms(C, n, csizes, 1, 10, 1, 0)
    print "A2\n"
    A2 = createA_onrouteoffroute_nonuniformcomms(C, n, csizes, 0.2, 5, 0.5, 0)
    exit()
    
    csizes = (4,4,4,4)
    n = sum(csizes)
    C = createC(n, csizes)
    print "A1\n"
    A1 = createA_fullcorridor_nonuniformcomms(C, n, csizes, 0.2, 5, 0.5, 1)
    print "A2\n"
    A2 = createA_fullcorridor_nonuniformcomms(C, n, csizes, 0.2, 5, 0.5, 0)
    exit()


    csizes = (3,3,3)
    n = sum(csizes)
    C = createC(n, csizes)
    A1, A2 = createA_simplecorridor_nonuniformcomms(C, n, csizes, 0.2, 100, 1, 0, 1)
    exit()

#    csizes = (5,5)
#    n = sum(csizes)
#    createA_weaklyconnectedcomms(n, csizes, 2, 0.8, 0.01, 10)
#    exit()

#    createA_simplecorridor_uniform(6, (2,2,2), 0.25, 0, 0.25)
    n = 9
    csizes = (n/3,n/3,n/3)
    C = createC(n, csizes)
    sdfactor = 0.01
    kfactor = 1000
    stage1CCIcrossk = 0.9
    stage2CCIcrossk = 0.0
    MIDcrossk = 0.3
    createA_simplecorridor_nonuniformK(C, n, csizes, sdfactor, kfactor, stage1CCIcrossk, stage2CCIcrossk, MIDcrossk)
    exit()
    
    n = 7
    createA_adopters_continuous(n, adoptercatgrowth=10.0)
    exit()

    csizes = (3,2)
    n = sum(csizes)
    C = createC(n, csizes)
    A = createA_uniformcomms(n, csizes, (0.2,0.2))

    csizes = (7,)
    n = sum(csizes)
    C = createC(n, csizes)
    A = createA_nonuniformK(C, n, csizes, sdfactor=0.01, kfactor=2, crosskfactor=10)
    print A
    print np.sum(A,axis=0)
    print np.sum(A,axis=1)

    csizes = (3,3)
    n = sum(csizes)
    C = createC(n, csizes)
    A = createA_uniformcomms(n, csizes, (0.3,0.3))
    print A

    csizes = (3,2)
    n = sum(csizes)
    C = createC(n, csizes)
    A = createA_uniformcomms(n, csizes, (0.1,0.1))
#    A = createA_weakandstrong(n, csizes, (0.3,0.5),0.2)
    print A
    for i in range(0,n):
        print np.sum(A[:,i])
    #print createApi(A,C)

    csizes = (6,4)
    n = sum(csizes)
    C = createC(n, csizes)
    A = createA_nonuniformcomms(C, n, csizes, 0.2, 500, 0.5)
    print A

if __name__ == "__main__":
    main()
