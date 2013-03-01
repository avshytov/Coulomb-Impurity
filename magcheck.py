import math
import numpy as np
import matplotlib
from matplotlib import pylab
from pylab import *
import diracsolver
from diracsolver import *
import ldos
from ldos import *

N = 300
rmin = 0.01
rmax = 25.0
r = zeros((N))
pot = zeros((N))
a = 0
mlist = zeros((2*a + 1))                                                        
#mlist = np.array(range(0,a))
mlist[0] = 0                                                                    
for i in range (0,N):
    r[i] = rmin +  i*(rmax-rmin) / N
#   pot[i] = -1.0 / 4.0 / r[i]
np.save("mlist",mlist)
np.save("potvec",pot)    

for B in [0.25, 0.50, 0.75, 1.0]:
    print "B=", B
    Emat, cdtens = diracham(r,pot,mlist,B)
    E, dostens = DOS(Emat, mlist ,r)

    ldosmat = ldoscalc()
    x = E
    lim = 3.3
    samps = [0.25, 0.5, 0.75, 1.0]
    for l in range (0, len(x)):
        if x[l-1] < -lim and x[l] >= - lim:
            lowl = l
        if x[l-1] < lim and x[l] >= lim:
            highl = l
            
    for rs in samps:
        si = int (rs / r[(len(r)-1)] * len(r))
        ri = r[si]
        y = ldosmat[:, si]
        print "r = ", r[si]
        for a in range (lowl, highl):
            if y[a-1] <= y[a] and y[a+1] < y[a]:
                print "peak position:", x[a], ", peak height:", y[a]
        plot(x,y, label='ldos slice as r = %g' %ri)
        xlim(x[lowl], x[highl])
legend()
show()
