from generateA import *
from calcTfromextensions import *
from calcE import *
from applyT import *
from updateG import *
import V2_3params
from threading import Thread
import numpy.linalg as LA

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
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

        #Add intra-cluster edges uniformly
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


#    pre_b_n = CHICAGO+JOLIET+DWIGHT+PONTIAC+CHENOA
#    plt.imshow(B_Nmat > 0, interpolation="nearest")

#    pre_kan = sum(OnRSIZES)+sum(OffRWSIZES)
#    kanmat = A[pre_kan:pre_kan+KANKAKEE,pre_kan:pre_kan+KANKAKEE]
#    axes = plt.gca()
#    axes.set_ylim([-0.5,kanmat.shape[0]-0.5])
#    axes.set_xlim([-0.5,kanmat.shape[1]-0.5])
#    plt.imshow(np.flip(np.sqrt(kanmat), axis=1), interpolation="nearest", cmap="viridis")
#    plt.title("Kankakee Adjacency Matrix")
#    plt.savefig("../plots/adjacency.png")
#    plt.show()

#    plt.imshow(B_N*100, interpolation="nearest")
#    plt.show()
#    plt.imshow(A, interpolation="nearest")
#    plt.show()

#NORMALIZE
    A += np.identity(sum(sizes))
    colsums = np.sum(A,axis=0)
    A/= colsums


    return A, sizes, NUMOnR



def iterate_full(C, commlangdistros, n, commsizes, alpha, Gfunc, iterations, Ts, presprob, G1 = None, A = None, EP = False, movements=None):
    cpoints = []
    cpointsnew = []
    for i in range(0,len(commsizes)):
        cpoints.append([])
        cpointsnew.append([])
    if G1 is None:
        G1 = Gfunc(C, commlangdistros) #initial lang distribution
        baseindex = 0
        for i, size in enumerate(commsizes):
            cpoints[i].append(np.mean(commlangdistros[baseindex:baseindex+size,:],axis=0))
            cpointsnew[i].append(np.mean(commlangdistros[baseindex:baseindex+size,:],axis=0))
            baseindex += size

    I = np.identity(A[:,0].size)
    invA = LA.inv(I-(1-alpha)*A)
    ages = np.ones(G1.shape[0])
    for i in range(0,iterations):
        G0 = G1
        if EP:
            E = calcE_geometric_EP(G0,C,A,alpha)
        else:
            E = calcE_geometric_noC_invA(G0,invA,alpha)
        G1 = applyTinf(E,Ts)
        G1, updatedindices = updateG_movement(C, G0, Gfunc, G1, presprob, movements) #update step
        ages *= (1-updatedindices)
        ages += (1-updatedindices)

        print ages[300:315]

        baseindex = 0
        for c, size in enumerate(commsizes):
            cpoints[c].append(np.mean(G1[baseindex:baseindex+size,:],axis=0))
            subG = G1[baseindex:baseindex+size,:]
            subyoung = ages[baseindex:baseindex+size]
            cpointsnew[c].append(np.mean(subG[subyoung > min(i-1,3)],axis=0))
            baseindex += size

    return cpoints, cpointsnew, G1



def iterate_learnthresh(C, commlangdistros, n, commsizes, alpha, Gfunc, iterations, p0, presprob, G1 = None, A = None, EP = False, movements=None):
    cpoints = []
    cpointsnew = []
    for i in range(0,len(commsizes)):
        cpoints.append([])
        cpointsnew.append([])
    if G1 is None:
        G1 = Gfunc(C, commlangdistros) #initial lang distribution
        baseindex = 0
        for i, size in enumerate(commsizes):
            cpoints[i].append(np.mean(commlangdistros[baseindex:baseindex+size,:],axis=0))
            cpointsnew[i].append(np.mean(commlangdistros[baseindex:baseindex+size,:],axis=0))
            baseindex += size

    I = np.identity(A[:,0].size)
    invA = LA.inv(I-(1-alpha)*A)
    ages = np.ones(G1.shape[0])
    for i in range(0,iterations):
        G0 = G1
        if EP:
            E = calcE_geometric_EP(G0,C,A,alpha)
        else:
            E = calcE_geometric_noC_invA(G0,invA,alpha)
        G1 = apply2varvariational_categorical(E,p0)
        G1, updatedindices = updateG_movement(C, G0, Gfunc, G1, presprob, movements) #update step
        ages *= (1-updatedindices)
        ages += (1-updatedindices)

#        print updatedindices.shape
        baseindex = 0
        for c, size in enumerate(commsizes):
            cpoints[c].append(np.mean(G1[baseindex:baseindex+size,:],axis=0))
            subG = G1[baseindex:baseindex+size,:]
            subyoung = ages[baseindex:baseindex+size]
            cpointsnew[c].append(np.mean(subG[subyoung > min(i-1,3)],axis=0))
            baseindex += size

    return cpoints, cpointsnew, G1




def plot_full_stlvs(allpoints, outname, STL, SPR, othercities, stages):
    for n, allpoint in enumerate(allpoints):
        outf = outname + ".png"
        print "Saving to", outf
#        for i, comm in enumerate(allpoint):
#            if i != STL and i != SPR and i not in othercities:
#                points = [point[0,0]*100 for point in comm]
#                lines = plt.plot(points)
#                plt.setp(lines, color='mediumpurple')
#                continue
#        for i, comm in enumerate(allpoint):
#            points = [point[0,0]*100 for point in comm]
#            if i == STL:
#                lines = plt.plot(points)
#                plt.setp(lines, color='red', linewidth = 3)
#            elif i == SPR:
#                lines = plt.plot(points)
#                plt.setp(lines, color='blue', linewidth = 3)
#            elif i in othercities:
#                lines = plt.plot(points)
#                plt.setp(lines, color='lightgreen')

        #STL
        points = [point[0,0]*100 for point in allpoint[STL]]
        lines = plt.plot(points)
        plt.setp(lines, color='red', linewidth = 3)
        #SPR
        points = [point[0,0]*100 for point in allpoint[SPR]]
        lines = plt.plot(points)
        plt.setp(lines, color='blue', linewidth = 3)

        axes = plt.gca()
        axes.set_ylim([0,100])
        axes.set_xlim([0,sum(stages)])
        
        #Plot other city and country averages
        sumpointsother = [0 for point in allpoint[1]]
        sumpointscountry = [0 for point in allpoint[1]]
        for i, comm in enumerate(allpoint):
            sumpointsother = [sumpointsother[j] + point[0,0]*100 if i in othercities else sumpointsother[j] for j, point in enumerate(comm)]
            avgpointsother = [sumpoint/len(othercities) for sumpoint in sumpointsother]
        avglineother = plt.plot(avgpointsother)
        plt.setp(avglineother, color='green', linewidth = 3)
        for i, comm in enumerate(allpoint):
            sumpointscountry = [sumpointscountry[j] + point[0,0]*100 if i not in othercities and i!=0 and i != len(allpoint)-1 and i!=STL and i!=SPR else sumpointscountry[j] for j, point in enumerate(comm)]
            avgpointscountry = [sumpoint/(len(allpoint)-len(othercities)-2) for sumpoint in sumpointscountry]
        avglinecountry = plt.plot(avgpointscountry)
        plt.setp(avglinecountry, color='purple', linewidth = 3)

        for i in range(len(stages)-1):
            plt.plot((sum(stages[0:i+1]),sum(stages[0:i+1])), (0, 100), 'k--', linewidth = 3)
        plt.title("Percent NCS+")
        plt.savefig(outf)
#        plt.show()
        plt.cla()


def plot_full_corridor(allpoints, outname, groupss, stages):

    def get_label(c, l):
        return mpatches.Patch(color=c, label=l)

    for n, allpoint in enumerate(allpoints):
        groups = groupss[n]
        leng0 = groups.count(0)-1
        leng1 = groups.count(1)-1
        outf = outname + ".png"
        print "Saving to", outf
        for i, comm in enumerate(allpoint):
#            if i == 0 or i == len(allpoint)-1:
#                continue
            points = [point[0,0]*100 for point in comm]
            lines = plt.plot(points)
            if i <= leng0:
                plt.setp(lines, color='salmon')
            else:
                plt.setp(lines, color='deepskyblue')
        axes = plt.gca()
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        axes.set_ylim([0,100])
        axes.set_xlim([0,sum(stages)])

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
        for i in range(len(stages)-1):
            plt.plot((sum(stages[0:i+1]),sum(stages[0:i+1])), (0, 100), 'k--', linewidth = 3)

        plt.ylabel('Percent NCS', fontsize = 15)
        plt.xlabel('Iteration #', fontsize = 15)
        plt.legend(handles=[get_label("red","On-Route Overall Avgs"),get_label("blue","Off-Route Overall Avgs"),get_label("salmon","On-Route Comm Avgs"),get_label("deepskyblue","Off-Route Comm Avgs"),mlines.Line2D([], [], linewidth=3, linestyle='--', color="black", label='Phase Boundary')], fontsize=10, loc="best")
        plt.title("Percent NCS across Iterations", fontsize=15)
        plt.tight_layout()
        plt.savefig(outf)
#        plt.show()
        plt.cla()


def corridor_geographic_migmigneut():
    stage1iters = 5
    stage2iters = 15
    stage3iters = 10
    stage4iters = 15
    stage5iters = 20
    presprob=0.9
    alpha = 0.2
    allpoints = []
    allpointsnew = []
    
    groupss = []
    moveg2 = np.matrix([[1.0,0.0]])
    moveg3 = np.matrix([[0.0,1.0]])
    moveg4 = np.matrix([[0.2,0.8]])


    print "St. Louis Corridor Neutral Change, Multi-Grammar Speakers"
    sdfactor = 0.2
    kfactor = 0.1
    intracrosskfactor = 0.05
    crosskfactor = 0.005
    Ts = [createT_advantage_alpha(0.3), createT_advantage_beta(0.3)]

    print "Creating Matrix 1"
    A1, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor, intracrosskfactor, True)
    print "Creating Matrix 2"
    A2, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor/10, intracrosskfactor, True)
    print "Creating Matrix 3"
    A3, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor, intracrosskfactor, False)
    print "Creating Matrix 4"
    A4, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor*3, intracrosskfactor, False)
    print "Creating Matrix 5"
    A5, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor*6, intracrosskfactor, False)
    n = sum(sizes)
    C = createC(n, (1,)*n)
    commlangdistros = np.repeat(np.matrix([[0.0, 1.0]]), n, axis=0)

    commlangdistros[0:sizes[0],:] = np.matrix([[1.0, 0.0]])

    depressionrate = 0.02
    postdepressionrate = 0.1
    midlandsrate = 0.1
    movements1 = [(0.1, moveg3) for i in xrange(sum(sizes))]
    movements2 = [(depressionrate, moveg2)  if i <= sum(sizes[0:NumOnR]) else  (midlandsrate, moveg3) for i in xrange(sum(sizes))]
    movements3 = [(postdepressionrate, moveg3) if i <= sum(sizes[0:NumOnR]) and i >= sizes[0] else (midlandsrate, moveg3) for i in xrange(sum(sizes))]
    movements4 = [(postdepressionrate, moveg3) if i <= sum(sizes[0:NumOnR]) and i >= sizes[0] else (depressionrate*3, moveg4) for i in xrange(sum(sizes))]
    movements5 = [(postdepressionrate, moveg3) if i <= sum(sizes[0:NumOnR]) and i >= sizes[0] else (midlandsrate, moveg3) for i in xrange(sum(sizes))]
    movechicago = ((0.0,moveg3),)*sizes[0]
    movements1[0:sizes[0]] = movechicago
    movements2[0:sizes[0]] = movechicago
    movements3[0:sizes[0]] = movechicago
    movements4[0:sizes[0]] = movechicago
    movements5[0:sizes[0]] = movechicago

    cpoints1, cpoints1new, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage1iters, Ts, presprob,  A=A1, movements=movements1)
    cpoints2, cpoints2new, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage2iters, Ts, presprob, G1=G, A=A2, movements=movements2)
    cpoints3, cpoints3new, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage3iters, Ts, presprob, G1=G, A=A3, movements=movements3)
    cpoints4, cpoints4new, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage4iters, Ts, presprob, G1=G, A=A4, movements=movements4)
    cpoints5, cpoints5new, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage5iters, Ts, presprob, G1=G, A=A5, movements=movements5)
    cpoints = cpoints1
    for i, points in enumerate(cpoints):
        points.extend(cpoints2[i])
    for i, points in enumerate(cpoints):
        points.extend(cpoints3[i])
    for i, points in enumerate(cpoints):
        points.extend(cpoints4[i])
    for i, points in enumerate(cpoints):
        points.extend(cpoints5[i])
    allpoints.append(cpoints)

    cpointsnew = cpoints1new
    for i, points in enumerate(cpointsnew):
        points.extend(cpoints2new[i])
    for i, points in enumerate(cpointsnew):
        points.extend(cpoints3new[i])
    for i, points in enumerate(cpointsnew):
        points.extend(cpoints4new[i])
    for i, points in enumerate(cpointsnew):
        points.extend(cpoints5new[i])
    allpointsnew.append(cpointsnew)


    groups = [0 if x+1 <= NumOnR else 1 for x in range(len(sizes))]
    groupss.append(groups)

    plot_full_corridor(allpoints, "../plots/stlgeo_migmig_neut", groupss, (stage1iters, stage2iters, stage3iters, stage4iters, stage5iters))
    plot_full_corridor(allpointsnew, "../plots/stlgeo_migmig_neut_youngonly", groupss, (stage1iters, stage2iters, stage3iters, stage4iters, stage5iters))
    plot_full_stlvs(allpoints, "../plots/stlgeo_migmig_neut_stlvs", STL=NumOnR-1, SPR=8, othercities=[5,3,15,18], stages=(stage1iters,stage2iters,stage3iters,stage4iters,stage5iters))


    
def corridor_geographic_migmigthresh(thresh):
    stage1iters = 5
    stage2iters = 15
    stage3iters = 10
    stage4iters = 15
    stage5iters = 20
    presprob=0.9
    alpha = 0.2
    allpoints = []
    allpointsnew = []
    
    groupss = []
    moveg2 = np.matrix([[1.0,0.0]])
    moveg3 = np.matrix([[0.0,1.0]])
    moveg4 = np.matrix([[0.2,0.8]])


    print "St. Louis Corridor Neutral Change, Multi-Grammar Speakers"
    sdfactor = 0.2
    kfactor = 0.1
    intracrosskfactor = 0.05
    crosskfactor = 0.01
    p0 = thresh
    print "Creating Matrix 1"
    A1, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor, intracrosskfactor, True)
    print "Creating Matrix 2"
    A2, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor/10, intracrosskfactor, True)
    print "Creating Matrix 3"
    A3, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor, intracrosskfactor, False)
    print "Creating Matrix 4"
    A4, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor*3, intracrosskfactor, False)
    print "Creating Matrix 5"
    A5, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor*6, intracrosskfactor, False)
    n = sum(sizes)
    C = createC(n, (1,)*n)
    commlangdistros = np.repeat(np.matrix([[0.0, 1.0]]), n, axis=0)

    commlangdistros[0:sizes[0],:] = np.matrix([[1.0, 0.0]])

    depressionrate = 0.02
    postdepressionrate = 0.0
    midlandsrate = 0.0
    movements1 = [(0.1, moveg3) for i in xrange(sum(sizes))]
    movements2 = [(depressionrate, moveg2)  if i <= sum(sizes[0:NumOnR]) else  (midlandsrate, moveg3) for i in xrange(sum(sizes))]
    movements3 = [(postdepressionrate, moveg3) if i <= sum(sizes[0:NumOnR]) and i >= sizes[0] else (midlandsrate, moveg3) for i in xrange(sum(sizes))]
    movements4 = [(postdepressionrate, moveg3) if i <= sum(sizes[0:NumOnR]) and i >= sizes[0] else (depressionrate*3, moveg4) for i in xrange(sum(sizes))]
    movements5 = [(postdepressionrate, moveg3) if i <= sum(sizes[0:NumOnR]) and i >= sizes[0] else (midlandsrate, moveg3) for i in xrange(sum(sizes))]
    movechicago = ((0.0,moveg3),)*sizes[0]
    movements1[0:sizes[0]] = movechicago
    movements2[0:sizes[0]] = movechicago
    movements3[0:sizes[0]] = movechicago
    movements4[0:sizes[0]] = movechicago
    movements5[0:sizes[0]] = movechicago

    print "\tMovement from CHI to On-Route", depressionrate

    cpoints1, cpoints1new, G = iterate_learnthresh(C, commlangdistros, n, sizes, alpha, initG_uniform, stage1iters, p0, presprob,  A=A1, movements=movements1)
    cpoints2, cpoints2new, G = iterate_learnthresh(C, commlangdistros, n, sizes, alpha, initG_uniform, stage2iters, p0, presprob, G1=G, A=A2, movements=movements2)
    print "\tMovement from Midlands to On-Route", postdepressionrate
    cpoints3, cpoints3new, G = iterate_learnthresh(C, commlangdistros, n, sizes, alpha, initG_uniform, stage3iters, p0, presprob, G1=G, A=A3, movements=movements3)
    cpoints4, cpoints4new, G = iterate_learnthresh(C, commlangdistros, n, sizes, alpha, initG_uniform, stage4iters, p0, presprob, G1=G, A=A4, movements=movements4)
    cpoints5, cpoints5new, G = iterate_learnthresh(C, commlangdistros, n, sizes, alpha, initG_uniform, stage5iters, p0, presprob, G1=G, A=A5, movements=movements5)
    cpoints = cpoints1
    for i, points in enumerate(cpoints):
        points.extend(cpoints2[i])
    for i, points in enumerate(cpoints):
        points.extend(cpoints3[i])
    for i, points in enumerate(cpoints):
        points.extend(cpoints4[i])
    for i, points in enumerate(cpoints):
        points.extend(cpoints5[i])
    allpoints.append(cpoints)

    cpointsnew = cpoints1new
    for i, points in enumerate(cpointsnew):
        points.extend(cpoints2new[i])
    for i, points in enumerate(cpointsnew):
        points.extend(cpoints3new[i])
    for i, points in enumerate(cpointsnew):
        points.extend(cpoints4new[i])
    for i, points in enumerate(cpointsnew):
        points.extend(cpoints5new[i])
    allpointsnew.append(cpointsnew)


    groups = [0 if x+1 <= NumOnR else 1 for x in range(len(sizes))]
    groupss.append(groups)

    plot_full_corridor(allpoints, "../plots/stlgeo_migmig_thresh_%s" % str(thresh).replace(".",""), groupss, (stage1iters, stage2iters, stage3iters, stage4iters, stage5iters))
    plot_full_corridor(allpointsnew, "../plots/stlgeo_migmig_thresh_%s_youngonly" % str(thresh).replace(".",""), groupss, (stage1iters, stage2iters, stage3iters, stage4iters, stage5iters))
    plot_full_stlvs(allpoints, "../plots/stlgeo_migmig_thresh_%s_stlvs" % str(thresh).replace(".",""), STL=NumOnR-1, SPR=8, othercities=[5,3,15,18], stages=(stage1iters,stage2iters,stage3iters,stage4iters,stage5iters))




def corridor_geographic_migdiffneut():
    stage1iters = 5
    stage2iters = 15
    stage3iters = 45
    stage4iters = 0
    stage5iters = 0
    presprob=0.9
    alpha = 0.2
    allpoints = []
    allpointsnew = []
    
    groupss = []
    moveg2 = np.matrix([[1.0,0.0]])
    moveg3 = np.matrix([[0.0,1.0]])
    moveg4 = np.matrix([[0.2,0.8]])


    print "St. Louis Corridor Neutral Change, Multi-Grammar Speakers"
    sdfactor = 0.2
    kfactor = 0.1
    intracrosskfactor = 0.02
    crosskfactor = 0.1
    Ts = [createT_advantage_alpha(0.3), createT_advantage_beta(0.3)]

    print "Creating Matrix 1"
    A1, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor, intracrosskfactor, True)
    print "Creating Matrix 2"
    A2, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor/10, intracrosskfactor, True)
    print "Creating Matrix 3"
    A3, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor, intracrosskfactor, False)
    n = sum(sizes)
    C = createC(n, (1,)*n)
    commlangdistros = np.repeat(np.matrix([[0.0, 1.0]]), n, axis=0)

    commlangdistros[0:sizes[0],:] = np.matrix([[1.0, 0.0]])

    depressionrate = 0.02
    postdepressionrate = 0.1
    midlandsrate = 0.01
    movements1 = [(0.1, moveg3) for i in xrange(sum(sizes))]
    movements2 = [(depressionrate, moveg2)  if i <= sum(sizes[0:NumOnR]) else  (midlandsrate, moveg3) for i in xrange(sum(sizes))]
    movements3 = [(postdepressionrate, moveg3) if i <= sum(sizes[0:NumOnR]) and i >= sizes[0] else (midlandsrate, moveg3) for i in xrange(sum(sizes))]
    movechicago = ((0.0,moveg3),)*sizes[0]
    movements1[0:sizes[0]] = movechicago
    movements2[0:sizes[0]] = movechicago
    movements3[0:sizes[0]] = movechicago

    cpoints1, cpoints1new, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage1iters, Ts, presprob,  A=A1, movements=movements1)
    cpoints2, cpoints2new, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage2iters, Ts, presprob, G1=G, A=A2, movements=movements2)
    cpoints3, cpoints3new, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage3iters, Ts, presprob, G1=G, A=A3, movements=movements3)
    cpoints = cpoints1
    for i, points in enumerate(cpoints):
        points.extend(cpoints2[i])
    for i, points in enumerate(cpoints):
        points.extend(cpoints3[i])
    allpoints.append(cpoints)

    cpointsnew = cpoints1new
    for i, points in enumerate(cpointsnew):
        points.extend(cpoints2new[i])
    for i, points in enumerate(cpointsnew):
        points.extend(cpoints3new[i])
    allpointsnew.append(cpointsnew)


    groups = [0 if x+1 <= NumOnR else 1 for x in range(len(sizes))]
    groupss.append(groups)

    plot_full_corridor(allpoints, "../plots/stlgeo_migdiff_neut", groupss, (stage1iters, stage2iters, stage3iters))
    plot_full_corridor(allpointsnew, "../plots/stlgeo_migdiff_neut_youngonly", groupss, (stage1iters, stage2iters, stage3iters))
    plot_full_stlvs(allpoints, "../plots/stlgeo_migdiff_neut_stlvs", STL=NumOnR-1, SPR=8, othercities=[5,3,15,18], stages=(stage1iters,stage2iters,stage3iters))







def corridor_geographic_migdiffadv():
    stage1iters = 5
    stage2iters = 15
    stage3iters = 20
    stage4iters = 25
    stage5iters = 0
    presprob=0.9
    alpha = 0.2
    allpoints = []
    allpointsnew = []
    
    groupss = []
    moveg2 = np.matrix([[1.0,0.0]])
    moveg3 = np.matrix([[0.0,1.0]])
    moveg4 = np.matrix([[0.2,0.8]])


    print "St. Louis Corridor Neutral Change, Multi-Grammar Speakers"
    sdfactor = 0.2
    kfactor = 0.1
    intracrosskfactor = 0.02
    crosskfactor = 0.1
    Ts = [createT_advantage_alpha(0.3), createT_advantage_beta(0.7)]

    print "Creating Matrix 1"
    A1, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor/10000, intracrosskfactor, True)
    print "Creating Matrix 2"
    A2, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor/10000, intracrosskfactor, True)
    print "Creating Matrix 3"
    A3, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor, intracrosskfactor, False)
    n = sum(sizes)
    C = createC(n, (1,)*n)
    commlangdistros = np.repeat(np.matrix([[0.0, 1.0]]), n, axis=0)

    commlangdistros[0:sizes[0],:] = np.matrix([[1.0, 0.0]])

    depressionrate = 0.02
    postdepressionrate = 0.1
    midlandsrate = 0.01
    movements1 = [(0.1, moveg3) for i in xrange(sum(sizes))]
    movements2 = [(depressionrate, moveg2)  if i <= sum(sizes[0:NumOnR]) else  (midlandsrate, moveg3) for i in xrange(sum(sizes))]
    movements3 = [(postdepressionrate, moveg3) if i <= sum(sizes[0:NumOnR]) and i >= sizes[0] else (midlandsrate, moveg3) for i in xrange(sum(sizes))]
    movements4 = [(postdepressionrate, moveg3) if i <= sum(sizes[0:NumOnR]) and i >= sizes[0] else (postdepressionrate, moveg3) for i in xrange(sum(sizes))]
    movechicago = ((0.0,moveg3),)*sizes[0]
    movements1[0:sizes[0]] = movechicago
    movements2[0:sizes[0]] = movechicago
    movements3[0:sizes[0]] = movechicago
    movements4[0:sizes[0]] = movechicago

    cpoints1, cpoints1new, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage1iters, Ts, presprob,  A=A1, movements=movements1)
    cpoints2, cpoints2new, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage2iters, Ts, presprob, G1=G, A=A2, movements=movements2)
    cpoints3, cpoints3new, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage3iters, Ts, presprob, G1=G, A=A3, movements=movements3)
    cpoints4, cpoints4new, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage4iters, Ts, presprob, G1=G, A=A3, movements=movements4)
    cpoints = cpoints1
    for i, points in enumerate(cpoints):
        points.extend(cpoints2[i])
    for i, points in enumerate(cpoints):
        points.extend(cpoints3[i])
    for i, points in enumerate(cpoints):
        points.extend(cpoints4[i])
    allpoints.append(cpoints)

    cpointsnew = cpoints1new
    for i, points in enumerate(cpointsnew):
        points.extend(cpoints2new[i])
    for i, points in enumerate(cpointsnew):
        points.extend(cpoints3new[i])
    for i, points in enumerate(cpointsnew):
        points.extend(cpoints4new[i])
    allpointsnew.append(cpointsnew)


    groups = [0 if x+1 <= NumOnR else 1 for x in range(len(sizes))]
    groupss.append(groups)

    plot_full_corridor(allpoints, "../plots/stlgeo_migdiff_adv", groupss, (stage1iters, stage2iters, stage3iters,stage4iters))
    plot_full_corridor(allpointsnew, "../plots/stlgeo_migdiff_adv_youngonly", groupss, (stage1iters, stage2iters, stage3iters,stage4iters))
    plot_full_stlvs(allpoints, "../plots/stlgeo_migdiff_adv_stlvs", STL=NumOnR-1, SPR=8, othercities=[5,3,15,18], stages=(stage1iters,stage2iters,stage3iters,stage4iters))
    



def corridor_geographic_diffneut():
    stage1iters = 5
    stage2iters = 15
    stage3iters = 45
    presprob=0.9
    alpha = 0.2
    allpoints = []
    allpointsnew = []
    
    groupss = []
    moveg2 = np.matrix([[1.0,0.0]])
    moveg3 = np.matrix([[0.0,1.0]])
    moveg4 = np.matrix([[0.2,0.8]])


    print "St. Louis Corridor Neutral Change, Multi-Grammar Speakers"
    sdfactor = 0.2
    kfactor = 0.1
    intracrosskfactor = 0.02
    crosskfactor = 0.1
    Ts = [createT_advantage_alpha(0.3), createT_advantage_beta(0.3)]

    print "Creating Matrix 1"
    A1, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor, intracrosskfactor, True)
    print "Creating Matrix 2"
    A2, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor/100, intracrosskfactor, True)
    print "Creating Matrix 3"
    A3, sizes, NumOnR  = createA_STLCorridor(sdfactor, kfactor, crosskfactor, intracrosskfactor, True)
    n = sum(sizes)
    C = createC(n, (1,)*n)
    commlangdistros = np.repeat(np.matrix([[0.0, 1.0]]), n, axis=0)

    commlangdistros[0:sizes[0],:] = np.matrix([[1.0, 0.0]])

    print "Running SImulation"
    depressionrate = 0.02
    postdepressionrate = 0.1
    midlandsrate = 0.01
    movements1 = [(0.0, moveg3) for i in xrange(sum(sizes))]
#    movements1 = [(0.0, moveg3) for i in xrange(sum(sizes))]
#    movements2 = [(depressionrate, moveg2)  if i <= sum(sizes[0:NumOnR]) else  (midlandsrate, moveg3) for i in xrange(sum(sizes))]
#    movements3 = [(postdepressionrate, moveg3) if i <= sum(sizes[0:NumOnR]) and i >= sizes[0] else (midlandsrate, moveg3) for i in xrange(sum(sizes))]
    movechicago = ((0.0,moveg3),)*sizes[0]
    movements1[0:sizes[0]] = movechicago
#    movements2[0:sizes[0]] = movechicago
#    movements3[0:sizes[0]] = movechicago

    cpoints1, cpoints1new, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage1iters, Ts, presprob,  A=A1, movements=movements1)
    cpoints2, cpoints2new, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage2iters, Ts, presprob, G1=G, A=A2, movements=movements1)
    cpoints3, cpoints3new, G = iterate_full(C, commlangdistros, n, sizes, alpha, initG_uniform, stage3iters, Ts, presprob, G1=G, A=A3, movements=movements1)
    cpoints = cpoints1
    for i, points in enumerate(cpoints):
        points.extend(cpoints2[i])
    for i, points in enumerate(cpoints):
        points.extend(cpoints3[i])
    allpoints.append(cpoints)

    cpointsnew = cpoints1new
    for i, points in enumerate(cpointsnew):
        points.extend(cpoints2new[i])
    for i, points in enumerate(cpointsnew):
        points.extend(cpoints3new[i])
    allpointsnew.append(cpointsnew)


    groups = [0 if x+1 <= NumOnR else 1 for x in range(len(sizes))]
    groupss.append(groups)

    plot_full_corridor(allpoints, "../plots/stlgeo_diff_neut", groupss, (stage1iters, stage2iters, stage3iters))
    plot_full_corridor(allpointsnew, "../plots/stlgeo_diff_neut_youngonly", groupss, (stage1iters, stage2iters, stage3iters))
    plot_full_stlvs(allpoints, "../plots/stlgeo_diff_neut_stlvs", STL=NumOnR-1, SPR=8, othercities=[5,3,15,18], stages=(stage1iters,stage2iters,stage3iters))


    
#dynamics_ncsadv_retreat()



    

corridor_geographic_diffneut()
corridor_geographic_migdiffneut()
corridor_geographic_migdiffadv()
corridor_geographic_migmigneut()
corridor_geographic_migmigthresh(0.3)
corridor_geographic_migmigthresh(0.5)
corridor_geographic_migmigthresh(0.8)



