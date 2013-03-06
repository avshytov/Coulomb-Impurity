import math
import numpy as np
import matplotlib
from matplotlib import pylab
from pylab import *
import diracsolver
from diracsolver import *
import ldos
from ldos import *

rmin = 0.01
rmax = 25.0
r = np.load("rvec.npy")
N = len(r)
pot = zeros((N))
a = 3
mlist = zeros((2*a + 1))
num = 7

mlist = np.array(range(0,a))
#mlist[0] = 0
for i in range (0,N):
   pot[i] = -0.9 / r[i]
np.save("mlist",mlist)
np.save("potvec",pot)    
Bs = np.arange(0.03,3.0,0.05)
samps = [1.0]# [0.25, 0.5, 0.75, 1.0]
lines = 11
allvals = np.zeros((lines,len(Bs), len(samps)))

for b in range (0, len(Bs)):
    B = Bs[b]
    print "B=", B
    Emat, cdtens = diracham(r,pot,mlist,B)
    E, dostens = DOS(Emat, mlist ,r)
    x = E
    ldosmat = ldoscalc()
    lim = 3.3
    for l in range (0, len(x)):
        if x[l-1] < -lim and x[l] >= - lim:
            lowl = l
        if x[l-1] < lim and x[l] >= lim:
            highl = l
    for c in range (0,len(samps)):    
        rs = samps[c]
        peaks = np.zeros((0))
        peaks = list(peaks)
        si = int (rs / r[(len(r)-1)] * len(r))
        ri = r[si]
        y = ldosmat[:, si]
        print "r = ", r[si]
        for a in range (lowl, highl):
            if y[a-1] <= y[a] and y[a+1] < y[a]:
                #print "peak position:", x[a], ", peak height:", y[a]
                peaks.extend([x[a]])
        np.save('peaks-B=%g-r=%g' %(B,ri), peaks)
        top = len(peaks) / 2 + lines / 2 + 1
        bot = len(peaks) / 2 - lines / 2
        allvals[:,b,c] = peaks[bot:top:]
np.save('peakpositions',allvals)
print allvals
for c in range (0, len(samps)):
    for a in range (0,lines):
        plot(Bs, allvals[a,:,c], label='%d th peak' %a)
        title('Energies vs. Field: r = %g')
        #figure()
        #plot(Bs, allvals[a,:,c]**2, label='%d th peak' %a)
        #title('Energies^2 vs. Field: r = %g' %r[samps[c]])
    show()

        
