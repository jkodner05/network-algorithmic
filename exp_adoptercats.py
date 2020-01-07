from generateA import *
from calcTfromextensions import *
from calcE import *
from applyT import *
from updateG import *
import V2_3params

import matplotlib.pyplot as plt


def iterate(C, commlangdistros, A, regenA, Gfunc, innovprob, Ts, alpha, presprob, iterations, catgrowth, C0, n, commsizes):
    G1 = Gfunc(C, commlangdistros, innovprob) #initial lang distribution
    points = [np.mean(G1,axis=0)]
    for i in range(0,iterations):
        G0 = G1
        if regenA:
            A = createA_adopters_continuous(n, catgrowth*n)
        E = calcE_geometric(G0,C,A,alpha)
        G1 = applyTneutral(E,Ts)
        G1 = updateG_preservepop(C, G0, Gfunc, G1, presprob, innovprob) #update step
        points.append(np.mean(G1,axis=0))
    return points

def plot(allpoints, outname, Ts):
    for plotg in range(0,len(Ts)):
        if plotg != 1:
            continue
        maximum = 0
        for i in range(0,len(allpoints)):
            points = [point[0,plotg] for point in allpoints[i]]
            maximum = max(maximum, max(points))
            lines = plt.plot(points)
            plt.setp(lines, color='b')
            axes = plt.gca()
            axes.set_ylim([0,1])
        plt.title("Grammar %s" % (plotg+1))
        if maximum > 0.0001:
            plt.savefig(outname + "_g%s.png" % (plotg+1))
            plt.show()
        else:
            plt.cla()

    plt.cla()
        
def vary_catgrowth():
    n = 200
    commsizes = (n,)
    plotg = 1
    iterations = 200
    C0 = createC(n, commsizes)
    C = createC(n, (1,)*n) #individual speaker level

    alpha = 0.90

    g = 2
    Ts = [createT_neutral(g) for i in range (0,g)]

    presprob = 0
    innovprob = 0.00
    allpoints = []
    for catgrowth in (0.0, 0.001, 0.002, 0.005, 0.01):#, 1.0, 10.0):#, 100.0, 1000.0):
        A = createA_adopters_continuous(n, catgrowth*n)
        commlangdistros = np.repeat(np.matrix([[1.0,0.0]]), n, axis=0)
        commlangdistros[0,:] = np.matrix([[0.0, 1.0]])
        print "Category Growth Factor: ", catgrowth
        points = iterate(C, commlangdistros, A, False, initG_uniform, innovprob, Ts, alpha, presprob, iterations, catgrowth, C0, n, commsizes)
        allpoints.append(points)
    plot(allpoints, "../plots/1adopters_catgrowth", Ts)


def vary_backadopters():
    n = 200
    commsizes = (n,)
    plotg = 1
    iterations = 200
    C0 = createC(n, commsizes)
    C = createC(n, (1,)*n) #individual speaker level

    alpha = 0.90

    g = 2
    Ts = [createT_neutral(g) for i in range (0,g)]

    presprob = 0
    innovprob = 0.00
    allpoints = []
    catgrowth = 0.01
    for catback in (0.0, 1, 5, 10, 25, 50):
        A = createA_adopters_continuous(n, catgrowth*n, catback)
        commlangdistros = np.repeat(np.matrix([[1.0,0.0]]), n, axis=0)
        commlangdistros[0,:] = np.matrix([[0.0, 1.0]])
        print "Gack Edge Factor", catback
        points = iterate(C, commlangdistros, A, False, initG_uniform, innovprob, Ts, alpha, presprob, iterations, catgrowth, C0, n, commsizes)
        allpoints.append(points)

    plot(allpoints, "../plots/1adopters_catback", Ts)



vary_catgrowth()
vary_backadopters()
