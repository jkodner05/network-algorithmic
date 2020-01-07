from generateA import *
from calcTfromextensions import *
from calcE import *
from applyT import *
from updateG import *
import V2_3params
import time

import matplotlib.pyplot as plt

alpha = 0.9
commlangdistros = np.matrix([[0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],[0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]])


def reset_seed():
    seed = "mustard"
    random.seed(seed)
    nprandom.seed(1234)



def iterate(C, commlangdistros, n, commsizes, crosscprob,alpha,Gfunc, iterations, Ts, presprob, A = None, EP = False):
    G1 = Gfunc(C, commlangdistros) #initial lang distribution
#    if A == None:
#        A = createA_uniformcomms(n, commsizes, (crosscprob,crosscprob))
    c1points = [commlangdistros[0,:]]
    c2points = [commlangdistros[1,:]]

    start = time.time()
    invApi = get_invApi(A, alpha)
    invC = get_invC(C)

    for i in range(0,iterations):
        G0 = G1
        if EP:
#            E = calcE_geometric_EP(G,C,A,alpha)
            E = calcE_geometric_EP_invApiinvC(G0,C,invC,invApi,alpha)
        else:
            E = calcE_geometric(G0,C,A,alpha)
        L = E.T
        G1 = updateG_preservepop(C, G0, Gfunc, L, presprob) #update step
        c1points.append(L[0,:])
        c2points.append(L[1,:])
    end = time.time()
    print end-start
    return c1points, c2points


def iterate_single(C, commlangdistros, n, commsizes, crosscprob,alpha,Gfunc, iterations, Ts, presprob, A = None, EP = False):
    G1 = Gfunc(C, commlangdistros) #initial lang distribution
#    if A == None:
#        A = createA_uniformcomms(n, commsizes, (crosscprob,crosscprob))
    c1points = [commlangdistros[0,:]]

    start = time.time()
    invApi = get_invApi(A, alpha)
    invC = get_invC(C)

    for i in range(0,iterations):
        G0 = G1
        if EP:
            #E = calcE_geometric_EP(G,C,A,alpha)
            E = calcE_geometric_EP_invApiinvC(G0,C,invC,invApi,alpha)
        else:
            E = calcE_geometric(G0,C,A,alpha)
        L = applyTinf(E,Ts)
        G1 = updateG_preservepop(C, G0, Gfunc, L, presprob) #update step
        c1points.append(L[0,:])
    end = time.time()
    print end-start

    return c1points


def plot_spec_gc(plotg, allpoints, outname):
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    maximum = 0
    for i in range(0,len(allpoints)):
        points = [point[0,plotg] for point in allpoints[i]]
        maximum = max(maximum, max(points))
        lines = plt.plot(points)
        plt.setp(lines, color='b')
        axes = plt.gca()
        axes.set_ylim([0,1])
    plt.title("Grammar 2", fontsize=22)
    plt.ylabel('Percent g2', fontsize = 20)
    plt.xlabel('Iteration #', fontsize = 20)
    plt.tight_layout()
    if maximum > 0.0001:
        plt.savefig(outname + "_g%s.png" % (plotg+1))
        with open(outname + "_g%s.txt" % (plotg+1), "w") as f:
            for i in range(0,len(allpoints)):
                points1 = [point[0,plotg] for point in allpoints[i]]
                f.write(" ".join([str(float(c)) for c in points1]))
                f.write("\n")
        plt.show()
    else:
        plt.cla()
    plt.cla()


def plot_spec_gc_withpred(plotg, allpoints, outname):
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    maximum = 0
    for i in range(0,len(allpoints)):
        points = [point[0,plotg] for point in allpoints[i]]
        maximum = max(maximum, max(points))
        lines = plt.plot(points)
        if i < len(allpoints)-1:
            plt.setp(lines, color='b')
        else:
            plt.setp(lines, color='r', linewidth = 4)
        axes = plt.gca()
        axes.set_ylim([0,1])

    plt.title("Grammar 2", fontsize=22)
    plt.ylabel('Percent g2', fontsize = 20)
    plt.xlabel('Iteration #', fontsize = 20)
    plt.tight_layout()
#    plt.title("Grammar %s" % (plotg))
    if maximum > 0.0001:
        plt.savefig(outname + "_g%s.png" % (plotg+1))
        with open(outname + "_g%s.txt" % (plotg+1), "w") as f:
            for i in range(0,len(allpoints)-1):
                points1 = [point[0,plotg] for point in allpoints[i]]
                f.write(" ".join([str(float(c)) for c in points1]))
                f.write("\n")
        plt.show()
    else:
        plt.cla()
    plt.cla()


def plot_histogram(plotg, allpoints, outname):
    maximum = 0
    points = [point[0,plotg-1] for point in allpoints]
    maximum = max(maximum, max(points))
    n, bins, patches = plt.hist(points, 25, facecolor="b", range=(0,1))
    print n
    print bins
    plt.title("Grammar %s @ 500th Generation" % (2))
    axes = plt.gca()
    if maximum > 0.0001:
        plt.savefig(outname + "_g%s.png" % (plotg))
        with open((outname + "_g%s.txt") % (plotg), "w") as f:
            for point in points:
                f.write(str(point) + "\n")
        with open((outname + "_binning_g%s.txt") % (plotg), "w") as f:
            f.write(" ".join([str(int(c)) for c in n]))
            f.write("\n")
            f.write(" ".join([str(float(c)) for c in bins]))
            f.write("\n")
    else:
        plt.cla()
    plt.cla()

    
def vary_Gupdate_EP():
    seed = "mustard"
    random.seed(seed)

    iterations = 500
    crosscprob = 0.3
    presprob = 0.0
    Ts = [createT(i, V2_3params.paramdict, V2_3params.sentencedict, m=10) for i in range (1, 9)]
    for n in (200, 2000, 20000):
        c1allpoints = []
        c2allpoints = []
        commsizes = (n/2,n/2)
        C = createC(n, commsizes)
        A = createA_uniformcomms(n, commsizes, (crosscprob,crosscprob))
        Api = createApi(A,C)
        for i in range(0, 10):
            print "\n\tEopulation Size: ", n
            c1points, c2points = iterate(C, commlangdistros, n, commsizes, crosscprob, alpha, initG_indicator_old, iterations, Ts, presprob, A=Api, EP=True)
            c1allpoints.append(c1points)
            c2allpoints.append(c2points)
        c1points = iterate_single(C, commlangdistros, n, commsizes, crosscprob, alpha, initG_uniform, iterations, Ts, presprob, A=Api, EP=True)
        c1allpoints.append(c1points)

        plot_spec_gc_withpred(5, c1allpoints, "../plots/2uniform_Gindicator_n%s_c1" % n)


def vary_Gupdate_EP_hist():
    seed = "mustard"
    random.seed(seed)

    iterations = 500
    for iterations in (1000,):
        crosscprob = 0.3
        presprob = 0.2
        Ts = [createT(i, V2_3params.paramdict, V2_3params.sentencedict, m=10) for i in range (1, 9)]
        for n in (200,2000,20000):
            commsizes = (n/2,n/2)
            C = createC(n, commsizes)
            A = createA_uniformcomms(n, commsizes, (crosscprob,crosscprob))
            Api = createApi(A,C)
            c1allpoints = []
            c2allpoints = []
            for i in range(0, 100):
                print "\n\tEopulation Size: ", n, "\tTrial: ", i
                c1points, c2points = iterate(C, commlangdistros, n, commsizes, crosscprob, alpha, initG_indicator_old, iterations, Ts, presprob, A=Api, EP=True)
                c1allpoints.append(c1points[-1])
                c2allpoints.append(c2points[-1])
            plot_histogram(5, c1allpoints, "../plots/2uniform_Gindicator_hist_%s_n%s_c1" % (iterations,n))


def vary_Gupdate_EP_advantage():
    seed = "mustard"
    random.seed(seed)

    commlangdistros = np.matrix([[0.98, 0.02]])
    iterations = 1000
    crosscprob = 0.5
    presprob = 0.2
    alpha = 0.31
    beta = 0.3
    Ts = [createT_advantage_alpha(alpha), createT_advantage_beta(beta)]
    for n in (200, 2000, 20000):
        c1allpoints = []
        commsizes = (n,)
        C = createC(n, commsizes)
        A = createA_uniformcomms(n, commsizes, (crosscprob,crosscprob))
        Api = createApi(A,C)
        i = 0
        while i < 10:
            print i
            print "\n\tEopulation Size: ", n
            c1points = iterate_single(C, commlangdistros, n, commsizes, crosscprob, alpha, initG_indicator_old, iterations, Ts, presprob, A=Api, EP=True)
            if c1points[iterations/20][0][:,1] == 0:
                print "redoing"
            else:
                i += 1
                c1allpoints.append(c1points)
        print "Eredicted"
        c1points = iterate_single(C, commlangdistros, n, commsizes, crosscprob, alpha, initG_uniform, iterations, Ts, presprob, A=Api, EP=True)
        c1allpoints.append(c1points)
        #plot_spec_gc(1, c1allpoints, "../plots/2uniform_Gindicator_adv_n%s_c1" % n)
        plot_spec_gc_withpred(1, c1allpoints, "../plots/2uniform_Gindicator_adv_n%s_c1" % n)


vary_Gupdate_EP()
reset_seed()
vary_Gupdate_EP_advantage()

#vary_Gupdate_EP_hist()
