import math
import numpy as np
import matplotlib
from matplotlib import pylab
from pylab import *

r = np.load("rvec.npy")
info = np.load("EMinfo.npy")
mlist = np.load("mlist.npy")
nm = len(mlist)
nr = len(r)
Ulist = [1e-3] #, 1e-5, 1e-4, 1e-3, 0.01, 0.1]
rlist = [r[16], r[32], r[48], r[64], r[80]]
drhoU = np.zeros((nr, nm, len(Ulist)))
tdrhoU = np.zeros((nr, len(Ulist)))
print Ulist

for a in range (0,len(Ulist)):
    print a
    drhoU[:,:,a] = np.load("mdrho-cou-U=%g-B=0-m=15-grid=400-E=-1--0.5.npy" %Ulist[a])
for m in range (0,nm):
        tdrhoU[:,:] += drhoU[:,m,:]  

#
# Plots dhro against r for all U in Ulist
#
figure()
for a in range (0,len(Ulist)):
    plot(r,tdrhoU[:,a], label="U = %g/r" %(-Ulist[a]))
    title("Change in charge density")
    legend()
show()

#
# Plots drho/Ustr against r for all U in Ulist
#
for a in range (0,len(Ulist)):
    plot(r,tdrhoU[:,a]/Ulist[a], label="U = %g/r" %(-Ulist[a]))
    title("Change in charge density divided by Ustr")
    legend()
show()

#
# Plots r resolved drho and plots against U
#
for b in range (0, len(rlist)):
    plot(Ulist, tdrhoU[b,:], '*', label="r = %g" %b)
    title("Change in Charge density against U for r =%g"%b)
    legend()
show()
