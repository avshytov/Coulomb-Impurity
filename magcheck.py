import math
import numpy as np
import matplotlib
from matplotlib import pylab
from pylab import *
import diracsolver
from diracsolver import *
import ldos
from ldos import *
import time

rmin = 0.01
rmax = 25.0
r = np.load("rvec-N=200.npy")
E = np.load("Evec-N=200.npy")
N = len(r)
pot = zeros((N))
info = np.zeros((2))
a = 1
Ustr = 1.5
info[0] = Ustr
#mlist = zeros((2*a + 1))
mlist = np.array(range(0,a))
mlist[0] = 0
for i in range (0,N):
   pot[i] = -Ustr / r[i]
np.save("mlist",mlist)
np.save("potvec",pot)    
Bs = [] 
for Q in range (0,106):
   Bs.extend([(0.03**2 * Q**2)])
np.save('Bs-%d' %len(Bs), Bs)
samps = [1.0]# [0.25, 0.5, 0.75, 1.0]
lines = 21
allvals = np.zeros((lines,len(Bs), len(samps)))
timestart = time.time()

for b in range (0, len(Bs)):
    timeB = time.time()
    print 'Current time elapsed:', (timeB - timestart)
    B = Bs[b]
    print "B=", B
    info[1] = B
    np.save("EMinfo", info)
    Ematp, Ematn, cdtens = diracham(r,pot,mlist,B)
    np.save('chargedenstens-N=%d-ms=%d-U=%g-B=%g' %(N,a,Ustr,B) ,cdtens)
    Es, dostens = DOS(Ematp, Ematn, mlist ,r)
    x = E
    ldosmat = ldoscalc()
#    ldosplot(ldosmat)
    lim = 21
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
       # np.save('peaks-B=%g-r=%g' %(B,ri), peaks)
        top = len(peaks) / 2 + lines / 2 + 1
        bot = len(peaks) / 2 - lines / 2
        allvals[:,b,c] = peaks[bot:top:]
np.save('peakpositions-Ustr=%g-%dms-grid=%d' %(Ustr, len(mlist), len(r)) ,allvals)
np.save('Evec.npy', E)
#print allvals
timeend = time.time()
print 'Total time taken:', (timeend - timestart)
for c in range (0, len(samps)):
    for a in range (0,lines):
        plot(Bs, allvals[a,:,c], '*', label='%d th peak' %a)
        title('Energies vs. Field: r = %g' %r[samps[c]])
figure()
for c in range (0, len(samps)):
    for a in range (0,lines):
        plot(Bs, allvals[a,:,c], label='%d th peak' %a)
        title('Energies vs. Field: r = %g' %r[samps[c]])
figure()
for c in range (0, len(samps)):
    for a in range (0,lines):
        plot(np.sqrt(Bs), allvals[a,:,c], '*', label='%d th peak' %a)
        title('Energies vs. sqrt(Field): r = %g' %r[samps[c]])
figure()
for c in range (0, len(samps)):
    for a in range (0,lines):
        plot(np.sqrt(Bs), allvals[a,:,c], label='%d th peak' %a)
        title('Energies vs. sqrt(Field): r = %g' %r[samps[c]])
legend()
show()

        
