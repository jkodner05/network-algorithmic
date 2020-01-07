from generateA import *
from calcTfromextensions import *
from calcE import *
from applyT import *
from updateG import *
import V2_3params

import matplotlib.pyplot as plt

import time

def reset_seed():
    seed = "mustard"
    random.seed(seed)
    nprandom.seed(1234)

def iterate_As(C, commlangdistros, n, commsizes, alpha, Gfunc, iterations, p0, presprob, As, G1=None):
    cpoints = []
    for i in range(0,len(commsizes)):
        cpoints.append([])
    if G1 is None:
        G1 = Gfunc(C, commlangdistros) #initial lang distribution
        baseindex = 0
        for i, size in enumerate(commsizes):
            cpoints[i].append(np.mean(commlangdistros[baseindex:baseindex+size,:],axis=0))
            baseindex += size

    print "Inverting As"
    start = time.time()
    invAs = [get_invA(A, alpha) for A in As]
    end = time.time()
    print end - start
    print "Running..."
    for i in range(0,iterations):
        print "iter", i
        start = time.time()
        G0 = G1

        E = calcE_geometric_noC_invA(G0,invAs[i % 10],alpha)

        G1 = apply2varvariational_categorical(E,p0)
        G1 = updateG_preservepop(C, G0, Gfunc, G1, presprob) #update step
        baseindex = 0
        for j, size in enumerate(commsizes):
            cpoints[j].append(np.mean(G1[baseindex:baseindex+size,:],axis=0))
            baseindex += size
        end = time.time()
        print end - start
    return cpoints, G1




def iterate_fullcc(C, commlangdistros, n, commsizes, alpha, Gfunc, iterations, p0, presprob, G1 = None, A = None):
    cpoints = []
    for i in range(0,len(commsizes)):
        cpoints.append([])
    if G1 is None:
        print "Creating G"
        G1 = Gfunc(C, commlangdistros) #initial lang distribution
        baseindex = 0
        for i, size in enumerate(commsizes):
            cpoints[i].append(np.mean(commlangdistros[baseindex:baseindex+size,:],axis=0))
            baseindex += size

    invA = get_invA(A, alpha)
    for i in range(0,iterations):
        print "iter", i
        start = time.time()
        G0 = G1

        E = calcE_geometric_noC_invA(G0,invA,alpha)
        G1 = apply2varvariational_categorical(E,p0)
        G1 = updateG_preservepop(C, G0, Gfunc, G1, presprob) #update step

        baseindex = 0
        for i, size in enumerate(commsizes):
            cpoints[i].append(np.mean(G1[baseindex:baseindex+size,:],axis=0))
            baseindex += size
        end = time.time()
        print end - start
    return cpoints, G1


def plot_full(allpoints, outname):
    axes = plt.gca()

    for n, allpoint in enumerate(allpoints):
        outf = outname + ("_%s.png" % (n*10))
        outfavg = outname + ("_%s_avg.png" % (n*10))
        print "Saving to", outf
        for i, comm in enumerate(allpoint):
            if i == 0 or i >= len(allpoint)-1:
                continue
            points = [100*point[0,0] for point in comm]
            lines = plt.plot(points)
        sumpoints = [point[0,0] for point in allpoint[1]]
        for i, comm in enumerate(allpoint):
            if i == 0 or i >= len(allpoint)-1 or i == 1:
                continue
            sumpoints = [sumpoints[j] + point[0,0] for j, point in enumerate(comm)]
        avgpoints = [100*sumpoint/(len(allpoint)-2) for sumpoint in sumpoints]


    plt.title("Percent Merged by Cluster", fontsize=22)
    tickstep = 20
    if len(avgpoints) <= 50:
        tickstep = 5
    elif len(avgpoints) <= 100:
        tickstep = 10
    plt.xticks(np.arange(0,len(avgpoints),step=tickstep),fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylabel('Percent Merged', fontsize = 20)
    plt.xlabel('Iteration #', fontsize = 20)

    axes.set_ylim([0,100])
    axes.set_xlim([0,len(avgpoints)-1])
    plt.tight_layout()
    plt.savefig(outf)
    avgline = plt.plot(avgpoints)
    plt.setp(avgline, color='black', linewidth = 4)
    plt.savefig(outfavg)
    plt.show()


def plot_avgs(allpoints, outname):
    for n, avgpoints in enumerate(allpoints):
        print outname
        outf = outname + ("_%s.png" % (n*10))
        print outf
        avgline = plt.plot([100*pt for pt in avgpoints], linewidth = 4)

    plt.title("Average Percent Merged", fontsize=22)
    tickstep = 20
    if len(avgpoints) <= 50:
        tickstep = 5
    elif len(avgpoints) <= 100:
        tickstep = 10
    plt.xticks(np.arange(0,len(avgpoints),step=tickstep),fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylabel('Percent Merged', fontsize = 20)
    plt.xlabel('Iteration #', fontsize = 20)
    axes = plt.gca()
    axes.set_ylim([0,100])
    axes.set_xlim([0,len(avgpoints)-1])
    plt.tight_layout()
    plt.savefig(outf)
    plt.show()


def merger_staticA(iters, numcomms, commsize):
    n = commsize*numcomms
    commsizes = (commsize,)*numcomms
    commlangdistros = np.repeat(np.matrix([[0.0, 1.0]]), n, axis=0)
    commlangdistros[commlangdistros[:,0].size-commsize:commlangdistros[:,0].size,:] = np.matrix([[1.0, 0.0]])
    C0 = createC(n, commsizes)
    C = createC(n, (1,)*n)
    presprob=0.8
    alpha = 0.45
    p0 = cotcaught_p0()
    allpoints = []

    A = createA_merger_nonuniformcomms(C0, n, commsizes, 0.2, 5, 0.5, 5)
    print "Running merger"
    cpoints, G = iterate_fullcc(C, commlangdistros, n, commsizes, alpha, initG_uniform, iters, p0, presprob,  A=A)
    allpoints.append(cpoints)

    plot_full(allpoints, "../plots/singleg_staticA%s_%s_%s" % (iters, numcomms, commsize))


def merger_staticA_nopres(iters, numcomms, commsize):
    n = commsize*numcomms
    commsizes = (commsize,)*numcomms
    commlangdistros = np.repeat(np.matrix([[0.0, 1.0]]), n, axis=0)
    commlangdistros[commlangdistros[:,0].size-commsize:commlangdistros[:,0].size,:] = np.matrix([[1.0, 0.0]])
    C0 = createC(n, commsizes)
    C = createC(n, (1,)*n)
    presprob=0.0
    alpha = 0.45
    p0 = cotcaught_p0()
    allpoints = []

    A = createA_merger_nonuniformcomms(C0, n, commsizes, 0.2, 3, 0.25, 5)
    print "Running merger"
    cpoints, G = iterate_fullcc(C, commlangdistros, n, commsizes, alpha, initG_uniform, iters, p0, presprob,  A=A)
    allpoints.append(cpoints)

    plot_full(allpoints, "../plots/singleg_staticA_nopres%s_%s_%s" % (iters, numcomms, commsize))


def merger_changingA(iters, numcomms, commsize):
    n = commsize*numcomms
    commsizes = (commsize,)*numcomms
    commlangdistros = np.repeat(np.matrix([[0.0, 1.0]]), n, axis=0)
    commlangdistros[commlangdistros[:,0].size-commsize:commlangdistros[:,0].size,:] = np.matrix([[1.0, 0.0]])
    C0 = createC(n, commsizes)
    C = createC(n, (1,)*n)
    presprob=0.9
    alpha = 0.45
    p0 = cotcaught_p0()
    allpoints = []

    As = []
    print "Preparing A matrices"
    for iterA in range(0,10):
        print iterA
        As.append(createA_merger_nonuniformcomms(C0, n, commsizes, 0.2, 5, 0.5, 5))

    print "Running merger"
    cpoints, G = iterate_As(C, commlangdistros, n, commsizes, alpha, initG_uniform, iters, p0, presprob,  As)
    allpoints.append(cpoints)

    plot_full(allpoints, "../plots/singleg_changingA%s_%s_%s" % (iters, numcomms, commsize))


def merger_avgs(iters, numcomms, commsize):
    n = commsize*numcomms
    commsizes = (commsize,)*numcomms
    commlangdistros = np.repeat(np.matrix([[0.0, 1.0]]), n, axis=0)
    commlangdistros[commlangdistros[:,0].size-commsize:commlangdistros[:,0].size,:] = np.matrix([[0.18, 0.82]])
    C0 = createC(n, commsizes)
    C = createC(n, (1,)*n)
    presprob=0.8
    alpha = 0.45
    allpoints = []

    trials = [cotcaught_p0()]*10
    for p0 in trials:
        print "Running iteration"
        print "Creating A"
        A1  = createA_merger_nonuniformcomms(C0, n, commsizes, 0.2, 5, 0.5, 5)
        print "Calculating Iteration"
        cpoints, G = iterate_fullcc(C, commlangdistros, n, commsizes, alpha, initG_uniform, iters, p0, presprob,  A=A1)

        sumpoints = [point[0,0] for point in cpoints[1]]
        for i, comm in enumerate(cpoints):
            print i
            if i == 0 or i == len(cpoints)-1 or i == 1:
                continue
            sumpoints = [sumpoints[j] + point[0,0] for j, point in enumerate(comm)]
        avgpoints = [sumpoint/(len(cpoints)-2) for sumpoint in sumpoints]

        allpoints.append(avgpoints)
        
    plot_avgs(allpoints, "../plots/singleg_avgs%s_%s_%s" % (iters, numcomms, commsize))


def merger_offset():
    iters = 100
    numcomms = 80
    commsize = 18
    n = commsize*numcomms
    commsizes = (commsize,)*numcomms
    commlangdistros = np.repeat(np.matrix([[0.0, 1.0]]), n, axis=0)
    commlangdistros[commlangdistros[:,0].size-commsize:commlangdistros[:,0].size,:] = np.matrix([[0.18, 0.82]])
    C0 = createC(n, commsizes)
    C = createC(n, (1,)*n)
    presprob=0.8
    alpha = 0.45
    allpoints = []



    fromCCI1 = 0
    fromCCI2 = 0
    p0 = cotcaught_p0()
    As = []
    print "Preparing A matrices"
    for iterA in range(0,iters):
        As.append(createA_fullcorridor_nonuniformcomms(C0, n, commsizes, 0.2, 5, 0.5, 0))

    for startoffset in (0,10,20):
        print "Running iteration, offset", startoffset
        cpoints, G = iterate_As(C, commlangdistros, n, commsizes, alpha, initG_uniform, iters, p0, presprob, As, startiteration=startoffset)

        sumpoints = [point[0,0] for point in cpoints[1]]
        for i, comm in enumerate(cpoints):
            if i == 0 or i == len(cpoints)-1 or i == 1:
                continue
            sumpoints = [sumpoints[j] + point[0,0] for j, point in enumerate(comm)]
        avgpoints = [sumpoint/(len(cpoints)-2) for sumpoint in sumpoints]

        allpoints.append(avgpoints)

    plot_avgs(allpoints, "../plots/singleg_offsets")

    
reset_seed()
merger_avgs(160,40,18)
exit()

reset_seed()
merger_staticA_nopres(20, 100, 75)

reset_seed()
merger_staticA(80, 40, 18)
merger_staticA(100, 100, 75)

reset_seed()
merger_changingA(120,100,75)
merger_changingA(120,40,18)


            
