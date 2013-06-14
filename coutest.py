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
Ulist = [-1.0, -0.1, 0.1, 1.0]
#Ulist = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4]
#Ulist = [-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1,
#          -0.05, -0.01, -0.001, 0.001, 0.01, 0.05,
#          0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
rlist = [16,32,48,64,80]
lim = len(Ulist) - 1
drhoU = np.zeros((nr, nm, len(Ulist)))
tdrhoU = np.zeros((nr, len(Ulist)))
linresp = np.load('generic-linear-response-cn-r0=1.0.npy')

for a in range (0,len(Ulist)):
    print "Loading drho for U=", Ulist[a]
    drhoU[:,:,a] = np.load("mdrho-cn-U=%g-B=0-m=15-grid=400-E=-1--0.5.npy" 
                           %Ulist[a])
for m in range (0,nm):
        tdrhoU[:,:] += drhoU[:,m,:]  
#
# Plots dhro against r for all U in Ulist
#
if True:
    for a in range (0,len(Ulist)/2):
        figure()
        plot(r,tdrhoU[:,a], label="U = %g/r" %(-Ulist[a]))
        plot(r,-tdrhoU[:,lim-a], label="U = %g/r 'flipped'" %(-Ulist[lim-a]))
        plot(r, Ulist[a] * linresp, label='Linear response')
        legend()
        ylabel("drho")
        xlabel("r")
        title("Change in charge density")
    show()
    for a in range (0,len(Ulist)/2):
        figure()
        loglog(r,tdrhoU[:,a], label="U = %g/r 'flipped'" %(-Ulist[a]))
        loglog(r,-tdrhoU[:,lim-a], label="U = %g/r" %(-Ulist[lim-a]))
        loglog(r, Ulist[a] * linresp, label='Linear response')
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
    for a in range (0,len(Ulist)/2):
        plot(r,tdrhoU[:,a]/Ulist[a], label="U = %g/r" %(-Ulist[a]))
        plot(r,tdrhoU[:,lim-a]/Ulist[lim-a], 
             label="U = %g/r" %(-Ulist[-(a+1)]))
    title("Change in charge density divided by Ustr")
    ylabel("drho/Ustr")
    xlabel("r")
    legend()
    show()
    figure()
    for a in range (0,len(Ulist)/2):
        loglog(r,-tdrhoU[:,a]/Ulist[a], 
               label="U = %g/r" %(-Ulist[a]))
        loglog(r,tdrhoU[:,a]/Ulist[lim-a], label="U = %g/r" %(-Ulist[lim-a]))
    loglog(r, r**(-1.0)/10.0, 'k--', label="Power -1")
    title("log plot of the change in charge density divided by Ustr")
    ylabel("drho/Ustr")
    xlabel("r")
    legend()
    show()
#
# Plots r resolved drho and plots against U
#
if True:
    Ulist = np.array(Ulist)
    figure()
    for b in rlist:
        plot(Ulist, tdrhoU[b,:], label="r = %g" %r[b])
    title("Change in Charge density against U for resolved r")
    ylabel("drho")
    xlabel("U")
    legend()
    figure()
    ra = 1e-3
    rb = 0.1
    for b in range (0,len(rlist)):
        ivals = [t[0] for t in enumerate(Ulist) if t[1] <= rb and t[1] >= ra]
        xvals = [Ulist[t] for t in ivals]
        yvals = [(abs(tdrhoU[b,t])) for t in ivals]
        grad = np.polyfit(log(xvals), log(yvals), 1)
        loglog(xvals,yvals, label="grad = %g" %grad[0])
        loglog(Ulist, -tdrhoU[b,:], label="r = %g" %r[b])    
        print "Gradient for r = ", r[rlist[b]], "Gradient =", grad[0]
#print "Gradient = ", grad[0], "for r=", r[rlist[b]]
    loglog(Ulist, Ulist/10.0, 'k--', label="linear")
    title("Absolute change in charge density against U for resolved r")
    ylabel("drho")
    xlabel("U")
    legend()
    show()

if True:
    Ulist = np.array(Ulist)
    print Ulist
    residue = np.zeros((nr,len(Ulist)))
    rval = [40,48]
    for a in range (0,len(Ulist)):
        residue[:,a] = tdrhoU[:,a] / Ulist[a] - linresp[:]
#    plot(r,residue[:,4], label='res')
#    plot(r, linresp, label='linresp')
#    plot(r,-tdrhoU[:,4] / Ulist[4], label='drho')
    legend()
    show()
    for b in rval:
        print abs(residue[b,:])
        plot(Ulist,residue[b,:], 'o-', label="r = %g" %r[b])
#        loglog(Ulist,residue[b,:], label="r = %g" %r[b])
    #loglog(Ulist,Ulist**2, '--', label="y = x^2")
    title("R resolved residue values against U")
    legend()
    show()
