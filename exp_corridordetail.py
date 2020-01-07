from generateA import *
from calcTfromextensions import *
from calcE import *
from applyT import *
from updateG import *
import V2_3params
from threading import Thread
import numpy.linalg as LA

import matplotlib.pyplot as plt
#np.set_printoptions(threshold=20)

def createA_STLCorridor(sdfactor, kfactor, crosskfactor, intracrosskfactor, isthegreatdepression):
    A = None
    sizefactor = 1.0 #1 = 1:200, 0.2 = 1:1000, 2 = 1:100, etc

    CHICAGO = 100
    MIDLANDS = 100

    JOLIET = 210
    DWIGHT = 10
    PONTIAC = 35
    CHENOA = 5
    B_N = 230
    ATLANTA = 5
    LEXINGTON = 75
    SPRINGFIELD = 380
    FARMERSVILLE = 5
    LITCHFIELD = 35
    MTOLIVE = 15
    COLLINSVILLE = 50
    STLOUIS = 4080
    OnRSIZES = [CHICAGO,JOLIET,DWIGHT,PONTIAC,CHENOA,B_N,ATLANTA,LEXINGTON,SPRINGFIELD,FARMERSVILLE,LITCHFIELD,MTOLIVE,COLLINSVILLE,STLOUIS]
    NUMOnR = len(OnRSIZES)

    OTTAWA = 80
    MINONK = 10
    ELPASO = 10
    PEORIA = 520
    HAVANA = 20
    PPLAINS = 5
    JACKSONVILLE = 100
    WHITEHALL = 15
    CARROLLTON = 10
    JERSEYVILLE = 25
    OffRWSIZES = [OTTAWA,MINONK,ELPASO,PEORIA,HAVANA,PPLAINS,JACKSONVILLE,WHITEHALL,CARROLLTON,JERSEYVILLE]
    NUMOffRW = len(OffRWSIZES)
    
    KANKAKEE = 110
    CLIFTON_ONARGA = 5
    PAXTON = 15
    U_C = 190
    TUSCOLA = 15
    ARGENTA = 5
    DECATUR = 295
    MATTOON = 80
    EFFINGHAM = 30
    VANDALIA = 25
    GREENVILLE = 15
    OffRESIZES = [KANKAKEE,CLIFTON_ONARGA,PAXTON,U_C,TUSCOLA,ARGENTA,DECATUR,MATTOON,EFFINGHAM,VANDALIA,GREENVILLE,MIDLANDS]
    NUMOffRE = len(OffRESIZES)

    sizes = copy(OnRSIZES)
    sizes.extend(OffRWSIZES)
    sizes.extend(OffRESIZES)
    sizes = [int(size*sizefactor) for size in sizes]

    connectpairs = []
    #OnR --> OnR connections
    for i in range(1,NUMOnR-1):
        connectpairs.append((i,i+1))
        connectpairs.append((i+1,i))
        if i+2 < NUMOnR:
            connectpairs.append((i,i+2))
            connectpairs.append((i+1,2))

    #Connect Chicago
    count = 0
    if isthegreatdepression:
        count = 10
    for i in range(count):
        #Chicago --> Joliet
        connectpairs.append((0,1))
        #Chicago --> Dwight
        connectpairs.append((0,2))
        #Chicago --> B_N
        connectpairs.append((0,5))
        #Chicago --> Springfield
        connectpairs.append((0,8))
        #Chicago --> STL
        connectpairs.append((0,NUMOnR-1))

    #Connect Midlands
    for i in range(1,len(sizes)-1):
        for j in range(10):
            connectpairs.append((len(sizes)-1, i))


    #OffRW --> OffRW connections
    for i in range(NUMOnR, NUMOnR+NUMOffRW-1):
        connectpairs.append((i,i+1))
        connectpairs.append((i+1,i))

    #OffRE --> OffRE connections
    for i in range(NUMOnR+NUMOffRW, NUMOnR+NUMOffRW+NUMOffRE-1):
        connectpairs.append((i,i+1))
        connectpairs.append((i+1,i))


    #OffRW --> OffRW connections
    for i in range(1, NUMOnR):
        if i + NUMOnR < len(sizes):
            connectpairs.append((i,i+NUMOnR))
            connectpairs.append((i+NUMOnR,i))
        if i + NUMOnR + NUMOffRW < len(sizes):
            connectpairs.append((i,i+NUMOnR+NUMOffRW))
            connectpairs.append((i+NUMOnR+NUMOffRW,i))


    #Joliet <--> Ottawa
    connectpairs.append((1,NUMOnR))
    connectpairs.append((NUMOnR,1))
    #Joliet <--> Kankakee
    connectpairs.append((1,NUMOnR+NUMOffRW))
    connectpairs.append((NUMOnR+NUMOffRW,1))

    #B_N <--> Peoria
    connectpairs.append((5,NUMOnR+3))
    connectpairs.append((NUMOnR+3,5))
    #B_N <--> U_C
    connectpairs.append((5,NUMOnR+NUMOffRW+3))
    connectpairs.append((NUMOnR+NUMOffRW+3,5))

    #Springfield <--> Jacksonville
    connectpairs.append((8,NUMOnR+6))
    connectpairs.append((NUMOnR+6,8))
    #Springfield <--> Decatur
    connectpairs.append((8,NUMOnR+NUMOffRW+6))
    connectpairs.append((NUMOnR+NUMOffRW+6,8))

    #Litchfield <--> Carrollton
    connectpairs.append((10,NUMOnR+8))
    connectpairs.append((NUMOnR+8,10))
    #Litchfield <--> Vandalia
    connectpairs.append((10,NUMOnR+NUMOffRW+9))
    connectpairs.append((NUMOnR+NUMOffRW+9,10))
    
    #STLouis <--> Jerseyville
    connectpairs.append((NUMOnR-1,NUMOnR+9))
    connectpairs.append((NUMOnR+9,NUMOnR-1))
    #STLouis <--> Greenville
    connectpairs.append((NUMOnR-1,NUMOnR+NUMOffRW+10))
    connectpairs.append((NUMOnR+NUMOffRW+10,NUMOnR-1))
    
    #for k,v in connectpairs:
    #    print k,v, "\t", sizes[k], sizes[v]

    n = sum(sizes)
    A = np.zeros((n,n))

    #Community internals
    #Divide up each community into clusters of <=20
    #Connections in each cluster distributed normally
    #Connect up the component clusters liberally
    subsizefactor = 20

    #Inter-cluster connections
    for fromc, toc in connectpairs:
        frombase = sum(sizes[:fromc])
        fromsize = sizes[fromc]
        tobase = sum(sizes[:toc])
        tosize = sizes[toc]
        Asub = A[frombase:frombase+fromsize, tobase:tobase+tosize]#np.zeros((fromsize,tosize))
        #Add inter-cluster edges uniformly
        k = int(((fromsize*tosize))*crosskfactor)

        iterfroms = [random.randint(0,(fromsize-1)/subsizefactor) for i in range(int(fromsize*crosskfactor)/5 + 1)]
        itertos = [random.randint(0,(tosize-1)/subsizefactor) for i in range(int(tosize*crosskfactor)/5 + 1)]
#        print len(iterfroms), len(itertos)
        for edge in range(k):
            iterfrom = random.choice(iterfroms)
            iterto = random.choice(itertos)
            i = -1
            while i >= min(fromsize,subsizefactor-1) or i < 0 or (20*iterfrom) + i >= fromsize:
                i = random.randint(0,min(fromsize,subsizefactor-1))
            j = -1
            while j >= min(tosize,subsizefactor-1) or j < 0 or (20*iterto) + j >= tosize:
                j = random.randint(0,min(tosize,subsizefactor-1))

            Asub[(20*iterfrom)+i,(20*iterto)+j] += 10000

        A[frombase:frombase+fromsize, tobase:tobase+tosize] = Asub

    baseindex = 0
    for c, csize in enumerate(sizes):
        Asub = np.zeros((csize,csize))

        #Add intra-cluster edges normally
        k = int(csize*csize*intracrosskfactor)
        mu = csize/2
        sd = csize*sdfactor/2
        for edge in range(k):
            i = -1
            while i >= csize or i < 0:
                i = int(round(random.gauss(mu, sd)))
#                i = random.randint(0,csize-1)
            j = i
            while j == i:
                j = random.randint(0,csize-1)
            Asub[i,j] += 100

        #divide cluster c of size c into subclusters size <= subsizefactor
        #overwrite intra-cluster edges within small clusters
        subsizes = [subsizefactor]*(csize/subsizefactor)
        if csize % subsizefactor:
            subsizes.append(csize % subsizefactor)
        subbaseindex = 0
        for s, subsize in enumerate(subsizes):
            Asubsub = np.zeros((subsize,subsize))

            #Create Gaussian cluster
            mu = subsize/2
            sd = subsize*sdfactor
            k = int(subsize*subsize*kfactor)
            for edge in range(k):
                i = -1
                while i >= subsize or i < 0:
                    i = int(round(random.gauss(mu, sd)))
                j = i
                while j == i:
                    j = random.randint(0,subsize-1)
                Asubsub[i,j] += 100

            Asub[subbaseindex:subbaseindex+subsize,subbaseindex:subbaseindex+subsize] = Asubsub
            subbaseindex += subsize
#        plt.imshow(Asub)
#        plt.show()


        
        A[baseindex:baseindex+csize,baseindex:baseindex+csize] = Asub
        baseindex += csize

        #NORMALIZE
    A += np.identity(sum(sizes))
    colsums = np.sum(A,axis=0)
    A/= colsums
#    plt.imshow(A > 0, cmap='gray')
#    plt.show()


    return A, sizes, NUMOnR



def iterate_full(C, commlangdistros, n, commsizes, alpha, Gfunc, iterations, Ts, presprob, G1 = None, A = None, EP = False, movements=None):
    cpoints = []
    for i in range(0,len(commsizes)):
        cpoints.append([])
    if G1 is None:
        G1 = Gfunc(C, commlangdistros) #initial lang distribution
        baseindex = 0
        for i, size in enumerate(commsizes):
            cpoints[i].append(np.mean(commlangdistros[baseindex:baseindex+size,:],axis=0))
            baseindex += size

    I = np.identity(A[:,0].size)
    invA = LA.inv(I-(1-alpha)*A)
    for i in range(0,iterations):
        G0 = G1
        if EP:
            E = calcE_geometric_EP(G0,C,A,alpha)
        else:
            E = calcE_geometric_noC_invA(G0,invA,alpha)
        G1 = applyTinf(E,Ts)
        if movements:
            G1 = updateG_movement(C, G0, Gfunc, G1, presprob, movements) #update step
        else:
            G1 = updateG_preservepop(C, G0, Gfunc, G1, presprob)

        baseindex = 0
        for i, size in enumerate(commsizes):
            cpoints[i].append(np.mean(G1[baseindex:baseindex+size,:],axis=0))
            baseindex += size
    return cpoints, G1


def iterate_fullcc(C, commlangdistros, n, commsizes, alpha, Gfunc, iterations, p0, presprob, G1 = None, A = None, EP = False, movements=None, neutral=False):
    cpoints = []
    for i in range(0,len(commsizes)):
        cpoints.append([])
    if G1 is None:
        G1 = Gfunc(C, commlangdistros) #initial lang distribution
        baseindex = 0
        for i, size in enumerate(commsizes):
            cpoints[i].append(np.mean(commlangdistros[baseindex:baseindex+size,:],axis=0))
            baseindex += size
    for i in range(0,iterations):
        G0 = G1
        if EP:
            E = calcE_geometric_EP(G0,C,A,alpha)
        else:
            E = calcE_geometric(G0,C,A,alpha)

        if not neutral:
            G1 = apply2varvariational_categorical(E,p0)
        else:
            G1 = applyTneutral(E,None)
        if movements:
            G1 = updateG_movement(C, G0, Gfunc, G1, presprob, movements) #update step
        else:
            G1 = updateG_preservepop(C, G0, Gfunc, G1, presprob)
        baseindex = 0
        for i, size in enumerate(commsizes):
            cpoints[i].append(np.mean(G1[baseindex:baseindex+size,:],axis=0))
            baseindex += size
    return cpoints, G


def plot_full_corridor(allpoints, outname, names, groupss, enddepression):
    for n, allpoint in enumerate(allpoints):
        groups = groupss[n]
        leng0 = groups.count(0)-1
        leng1 = groups.count(1)-1
        outf = outname + ("_%s.png" % names[n])
        print "Saving to", outf
        for i, comm in enumerate(allpoint):
#            if i == 0 or i == len(allpoint)-1:
#                continue
            points = [point[0,0]*100 for point in comm]
            lines = plt.plot(points)
            if i <= leng0:
                plt.setp(lines, color='magenta')
            else:
                plt.setp(lines, color='cyan')
            axes = plt.gca()
            axes.set_ylim([0,100])
        sumpoints = [0 for point in allpoint[1]]
        sumpoints0 = [0 for point in allpoint[1]]
        sumpoints1 = [0 for point in allpoint[1]]
        for i, comm in enumerate(allpoint):
            if i == 0 or i == len(allpoint)-1:
                continue
            sumpoints = [sumpoints[j] + point[0,0]*100 for j, point in enumerate(comm)]
            sumpoints0 = [sumpoints0[j] + point[0,0]*100 if i <= leng0 else sumpoints0[j] for j, point in enumerate(comm)]
            sumpoints1 = [sumpoints1[j] + point[0,0]*100 if i > leng0 else sumpoints1[j] for j, point in enumerate(comm)]
        avgpoints = [sumpoint/(len(allpoint)-2) for sumpoint in sumpoints]
        if leng0 > 0:
            avgpoints0 = [sumpoint/leng0 for sumpoint in sumpoints0]
        else:
            avgpoints0 = [0 for point in allpoint[1]]
        if leng1 > 0:
            avgpoints1 = [sumpoint/leng1 for sumpoint in sumpoints1]
        else:
            avgpoints1 = [0 for point in allpoint[1]]
        avgline0 = plt.plot(avgpoints0)
        avgline1 = plt.plot(avgpoints1)
        plt.setp(avgline0, color='red', linewidth = 3)
        plt.setp(avgline1, color='blue', linewidth = 3)
        depressionline = plt.plot((enddepression,enddepression), (0, 100), 'k--', linewidth = 3)
        plt.title("Percent NCS+")
        plt.savefig(outf)
        plt.show()


def plot_full_stlvs(allpoints, outname, names, STL, SPR, othercities, enddepression):
    for n, allpoint in enumerate(allpoints):
        outf = outname + ("_%s.png" % names[n])
        print "Saving to", outf
        for i, comm in enumerate(allpoint):
            points = [point[0,0]*100 for point in comm]
            if i == STL:
                lines = plt.plot(points)
                plt.setp(lines, color='red')
            elif i == SPR:
                lines = plt.plot(points)
                plt.setp(lines, color='blue')
            elif i in othercities:
                lines = plt.plot(points)
                plt.setp(lines, color='green')
            else:
                continue
            axes = plt.gca()
            axes.set_ylim([0,100])
        sumpointsother = [0 for point in allpoint[1]]
        for i, comm in enumerate(allpoint):
            sumpointsother = [sumpointsother[j] + point[0,0]*100 if i in othercities else sumpointsother[j] for j, point in enumerate(comm)]
            avgpointsother = [sumpoint/len(othercities) for sumpoint in sumpointsother]
        avglineother = plt.plot(avgpointsother)
        plt.setp(avglineother, color='green', linewidth = 3)
        depressionline = plt.plot((enddepression,enddepression), (0, 100), 'k--', linewidth = 3)
        plt.title("Percent NCS+")
        plt.savefig(outf)
        plt.show()


def dynamics_ncsadv_retreat():
    stage1iters = 5
    stage3iters = 10
    stage4iters = 10
    presprob=0.9
    alpha = 0.2
    allpoints = []
    advantages = ((0.82,0.82),(0.8,0.82),(0.76,0.82))
    advantages = ((0.82,0.82),)
    
    names = []
    groupss = []
    moveg2 = np.matrix([[1.0,0.0]])
    moveg3 = np.matrix([[0.0,1.0]])


    print "St. Louis Corridor Neutral Change, Multi-Grammar Speakers"
    sdfactor = 0.2
    kfactor = 0.1
    intracrosskfactor = 0.02
    crosskfactor = 0.1
    for a, b in advantages:
        Ts = [createT_advantage_alpha(a), createT_advantage_beta(b)]
        names.append(str(b-a))
        print "a", a
        print "b", b
        for stage2iters in (5,)*1:

#            print "Creating Matrix 1"
            A1, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor/10, intracrosskfactor, True)
            print "Creating Matrix 2"
            A2, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor/10, intracrosskfactor, True)
            print "Creating Matrix 3"
            A3, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor, intracrosskfactor, False)
            print "Creating Matrix 4"
            A4, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor*10, intracrosskfactor, False)
            print "Creating Matrix 5"
            A5, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor/10, intracrosskfactor, False)
            n = sum(sizes)
            C = createC(n, (1,)*n)
            commlangdistros = np.repeat(np.matrix([[0.0, 1.0]]), n, axis=0)

            commlangdistros[0:sizes[0],:] = np.matrix([[1.0, 0.0]])

            depressionrate = 0.05
            postdepressionrate = 0.1
            midlandsrate = 0.01
            #exclude CHI
            movements2 = [(depressionrate, moveg2)  if i <= sum(sizes[0:NumOnR]) else  (midlandsrate, moveg3) for i in xrange(sum(sizes))]
            #exclude CHI and STL
            movements3 = [(postdepressionrate, moveg3) if i <= sum(sizes[0:NumOnR])-5 and i >= sizes[0] else (midlandsrate, moveg3) for i in xrange(sum(sizes))]
            movements4 = [(postdepressionrate, moveg3) if i <= sum(sizes[0:NumOnR])-5 and i >= sizes[0] else (midlandsrate, moveg3) for i in xrange(sum(sizes))]
            movechicago = ((0.0,moveg3),)*sizes[0]
            movestl = ((0.05,moveg2),)*sizes[NumOnR-1]
            movements3[0:sizes[0]] = movechicago
            movements4[0:sizes[0]] = movechicago
            movements3[sum(sizes[:NumOnR-1]):sum(sizes[:NumOnR])] = movestl

            print "\tMovement from CHI to On-Route", depressionrate
            cpoints, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage1iters, Ts, presprob,  A=A1)
            cpoints2, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage2iters, Ts, presprob, G1=G, A=A2, movements=movements2)
#            cpoints2, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage2iters, Ts, presprob, A=A2, movements=movements2)
            print "\tMovement from Midlands to On-Route", postdepressionrate
            cpoints3, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage3iters, Ts, presprob, G1=G, A=A3, movements=movements3)
            cpoints4, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage4iters/2, Ts, presprob, G1=G, A=A4, movements=movements4)
            cpoints5, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage4iters/2, Ts, presprob, G1=G, A=A5, movements=movements4)
#            cpoints = cpoints2
            for i, points in enumerate(cpoints):
                points.extend(cpoints2[i])
            for i, points in enumerate(cpoints):
                points.extend(cpoints3[i])
            for i, points in enumerate(cpoints):
                points.extend(cpoints4[i])
            for i, points in enumerate(cpoints):
                points.extend(cpoints5[i])

            allpoints.append(cpoints)

            groups = [0 if x+1 <= NumOnR else 1 for x in range(len(sizes))]
            groupss.append(groups)

    plot_full_corridor(allpoints, "../plots/stldetail_ncsadv_retreat", names, groupss, stage1iters+stage2iters)
    plot_full_stlvs(allpoints, "../plots/stldetail_ncsadv_stlvs", names, STL=NumOnR-1, SPR=8, othercities=[5,3,15,18], enddepression=stage1iters+stage2iters)


    
dynamics_ncsadv_retreat()

