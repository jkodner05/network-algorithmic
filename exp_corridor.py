from generateA import *
from calcTfromextensions import *
from calcE import *
from applyT import *
from updateG import *
import V2_3params

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

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
#            E = calcE_geometric(G0,C,A,alpha)
        G1 = applyTinf(E,Ts)
        if movements:
            G1, updatedindices = updateG_movement(C, G0, Gfunc, G1, presprob, movements) #update step
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
            G1, updatedindices = updateG_movement(C, G0, Gfunc, G1, presprob, movements) #update step
        else:
            G1 = updateG_preservepop(C, G0, Gfunc, G1, presprob)
        baseindex = 0
        for i, size in enumerate(commsizes):
            cpoints[i].append(np.mean(G1[baseindex:baseindex+size,:],axis=0))
            baseindex += size
    return cpoints, G


def plot_full_corridor(allpoints, outname, names, groupss, enddepression):

    def get_label(c, l):
        return mpatches.Patch(color=c, label=l)

    for n, allpoint in enumerate(allpoints):
        groups = groupss[n]
        leng0 = groups.count(0)-1
        leng1 = groups.count(1)-1
        outf = outname + ("_%s.png" % names[n])
        print "Saving to", outf
        for i, comm in enumerate(allpoint):
            if i == 0 or i == len(allpoint)-1:
                continue
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
        axes.set_xlim([0,len(comm)-1])
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
        plt.ylabel('Percent NCS', fontsize = 15)
        plt.xlabel('Iteration #', fontsize = 15)

        depressionline = plt.plot((enddepression,enddepression), (0, 100), 'k--', linewidth = 3)
        plt.legend(handles=[get_label("red","On-Route Overall Avgs"),get_label("blue","Off-Route Overall Avgs"),get_label("salmon","On-Route Comm Avgs"),get_label("deepskyblue","Off-Route Comm Avgs"),mlines.Line2D([], [], linewidth=3, linestyle='--', color="black", label='Phase Boundary')], fontsize=10, loc="best")
        plt.title("Percent NCS across Iterations", fontsize=15)
        plt.tight_layout()
        plt.savefig(outf)
        plt.show()


def dynamics_ncsadv_retreat():
    stage1iters = 5
    stage3iters = 75
    numcomms = 40
    commsize = 18
    n = commsize*numcomms
    commsizes = (commsize,)*numcomms
    commlangdistros = np.repeat(np.matrix([[0.0, 1.0]]), n, axis=0)
    commlangdistros[0:commsize,:] = np.matrix([[1.0, 0.0]])

    C0 = createC(n, commsizes)
    C = createC(n, (1,)*n)

    presprob=0.9
    alpha = 0.45
    allpoints = []
    advantages = ((0.82,0.76),(0.82,0.8),(0.82,0.82),(0.8,0.82),(0.76,0.82))
    
    fromCHI1 = 0
    names = []
    groupss = []
    moveg2 = np.matrix([[1.0,0.0]])
    moveg3 = np.matrix([[0.0,1.0]])
    fromCHI2 = numcomms/2

    print "St. Louis Corridor Neutral Change, Multi-Grammar Speakers"
    for a, b in advantages:
        Ts = [createT_advantage_alpha(a), createT_advantage_beta(b)]
        names.append(str(b-a))
        print "a", a
        print "b", b
        for stage2iters in (20,)*1:

            A1  = createA_onrouteoffroute_nonuniformcomms(C0, n, commsizes, 0.2, 5, 0.5, fromCHI1, True)
            A2  = createA_onrouteoffroute_nonuniformcomms(C0, n, commsizes, 0.2, 5, 0.5, fromCHI2, True)
            A3  = createA_onrouteoffroute_nonuniformcomms(C0, n, commsizes, 0.2, 5, 0.5, fromCHI1, True)

            depressionrate = 0.02
            postdepressionrate = 0.05
            movements2 = [(depressionrate, moveg2)  if i <= fromCHI2*commsize else  (0.00, moveg3) for i in xrange(n)]
            movements3 = [(postdepressionrate, moveg3) if i <= fromCHI2*commsize and i > 0 else (0.00, moveg3) for i in xrange(n)]

            print "\tMovement from CHI to On-Route", depressionrate
            cpoints, G = iterate_full(C, commlangdistros, n, commsizes, alpha, initG_uniform, stage1iters, Ts, presprob,  A=A1)
            cpoints2, G = iterate_full(C, commlangdistros, n, commsizes, alpha, initG_uniform, stage2iters, Ts, presprob, G1=G, A=A2, movements=movements2)
            print "\tMovement from Midlands to On-Route", postdepressionrate
            cpoints3, G = iterate_full(C, commlangdistros, n, commsizes, alpha, initG_uniform, stage3iters, Ts, presprob, G1=G, A=A3, movements=movements3)
            for i, points in enumerate(cpoints):
                points.extend(cpoints2[i])
            for i, points in enumerate(cpoints):
                points.extend(cpoints3[i])
            allpoints.append(cpoints)

            groups = [0 if x+1 <= fromCHI2 else 1 for x in range(numcomms)]
            groupss.append(groups)

    plot_full_corridor(allpoints, "../plots/dynamics_ncsadv_retreat", names, groupss, stage1iters+stage2iters)


    
dynamics_ncsadv_retreat()
