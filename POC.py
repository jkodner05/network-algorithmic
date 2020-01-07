import numpy as np
import numpy.linalg as LA

#constants
H = np.matrix('1 0 ; 1 0 ; 1 0 ; 0 1 ; 0 1')
Hplus = LA.inv(np.transpose(H)*H)
B_0 = np.matrix('0 1 ; 1 0 ; 1 0 ; 0 1 ; 0 1')
n = 5
c = 2
g = 2
alpha = 0.1

#Admits EP
A_0 = np.matrix('.0 .4 .4 .1 .1 ; .4 .0 .4 .1 .1 ; .4 .4 .0 .1 .1 ; .1 .1 .1 .0 .7 ; .1 .1 .1 .7 .0')
T_0 = np.matrix('1 .3 ; 0 .7')
T_1 = np.matrix('.2 0 ; .8 1')

I = np.identity(n)
Ipi = np.identity(c)

#Look at this. It's interpretable 
Api = Hplus * np.transpose(H) * A_0 * H

print "Api:\n", Api, "\n"

#should be equal and column stochastic
Ppi_1 = alpha * np.transpose(B_0) * H * LA.inv(Ipi-(1-alpha)*Api) * Hplus
P_1 = np.transpose(B_0) * alpha * LA.inv(I-(1-alpha)*A_0) * H * Hplus

print "Ppi_1 == P_1:", np.allclose(Ppi_1, P_1)

print "Ppi_1 col stoch:", np.allclose(np.sum(Ppi_1, axis=0),np.matrix('1 1'))


w, e_0 = LA.eig(P_1[0,0] * T_0 + P_1[1,0] * T_1)
domeigv = np.where(abs(w-1) < 0.000001)[0][0]
L_10 = e_0[:,domeigv]/np.sum(e_0[:,domeigv])
print "L_10: Gen 1 lang distro in community 0"
print L_10

w, e_1 = LA.eig(P_1[0,1] * T_0 + P_1[1,1] * T_1)
domeigv = np.where(abs(w-1) < 0.000001)[0][0]
L_11 = e_1[:,domeigv]/np.sum(e_1[:,domeigv])
print "L_10: Gen 1 lang distro in community 1"
print L_11
