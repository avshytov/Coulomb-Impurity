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
Ulist = [1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 
         3e-3, 5e-3, 0.01, 0.03, 0.05, 0.075, 0.1]
rlist = [r[12], r[16], r[32], r[48], r[64], r[80]]
drhoU = np.zeros((nr, nm, len(Ulist)))
tdrhoU = np.zeros((nr, len(Ulist)))

for a in range (0,len(Ulist)):
    print "Loading drho for U=", Ulist[a]
    drhoU[:,:,a] = np.load("mdrho-cou-U=%g-B=0-m=15-grid=400-E=-1--0.5.npy" %Ulist[a])
for m in range (0,nm):
        tdrhoU[:,:] += drhoU[:,m,:]  

#
# Plots dhro against r for all U in Ulist
#
if True:
    figure()
    for a in range (0,len(Ulist)):
        plot(r,tdrhoU[:,a], label="U = %g/r" %(-Ulist[a]))
    ylabel("drho")
    xlabel("r")
    title("Change in charge density")
    legend()
    show()
    figure()
    for a in range (0,len(Ulist)):
        loglog(r,-tdrhoU[:,a], label="U = %g/r" %(-Ulist[a]))
    title("log plot of the change in charge density")
    ylabel("drho")
    xlabel("r")
    legend()
    show()
#
# Plots drho/Ustr against r for all U in Ulist
#
if True:
    figure()
    for a in range (0,len(Ulist)):
        plot(r,tdrhoU[:,a]/Ulist[a], label="U = %g/r" %(-Ulist[a]))
    title("Change in charge density divided by Ustr")
    ylabel("drho/Ustr")
    xlabel("r")
    legend()
    show()
    figure()
    for a in range (0,len(Ulist)):
        loglog(r,-tdrhoU[:,a]/Ulist[a], label="U = %g/r" %(-Ulist[a]))
    title("log plot of the change in charge density divided by Ustr")
    ylabel("drho/Ustr")
    xlabel("r")
    legend()
    show()
#
# Plots r resolved drho and plots against U
#
if True:
    figure()
    for b in range (0, len(rlist)):
        plot(Ulist, tdrhoU[b,:], label="r = %g" %rlist[b])
    title("Change in Charge density against U for r =%g"%rlist[b])
    ylabel("drho")
    xlabel("U")
    legend()
    figure()
    for b in range (0,len(rlist)):
        loglog(Ulist, -tdrhoU[b,:], label="r = %g" %rlist[b])
    loglog(Ulist, Ulist, 'k--', label="linear")
    title("Change in Charge density against U for r =%g"%rlist[b])
    ylabel("drho")
    xlabel("U")
    legend()
    show()
