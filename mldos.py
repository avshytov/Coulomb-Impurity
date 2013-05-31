import numpy as np
import math
import matplotlib
from matplotlib import pylab
from pylab import *

def mldoscalc():
    r = np.load("rvec.npy")
    info = np.load("EMinfo.npy")
    mlist = np.load("mlist.npy")
    nm = len(mlist)
    nr = len(r)
    dostens = np.load("dostens-U=%g-B=%g-ms=%d-N=%d.npy" %(info[0], info[1], nm, nr))
    cdtens = np.load("cdtens-U=%g-B=%g-ms=%d-N=%d.npy" %(info[0], info[1],nm,nr))
    ni = len(dostens[0,0])
    nn = len(dostens[0])
    mldostens = np.zeros((nm,ni,nr))
    cdmat = np.zeros((nr,nn))
    dosmat = np.zeros((nn,ni))
    ldosmat = np.zeros((nr,ni))
    for m in range (0,nm):
        print "Calculating mLDOS for:", (m+1), "/", nm
        cdmat = cdtens[m,:,:]
        dosmat = dostens[m,:,:]
        ldosmat = np.dot(cdmat, dosmat)
        ldosmat = np.transpose(ldosmat)
        mldostens[m,:,:] = ldosmat[:,:]
    print "Saving data: U=", info[0], "ms = ", nm, "B = ", info[1], "grid = ", nr
    print "max, min: ", np.max(mldostens), np.min(mldostens)
    np.save("mldostens-U=%g-ms=%d-B=%g-grid=%d" %(info[0], nm, info[1], nr), mldostens)
    return mldostens

def mcdplotter(mldostens):
    r = np.load("rvec.npy")
    mlist = np.load("mlist.npy")
    E = np.load("Evec.npy")
    info = np.load("EMinfo.npy")
    nr = len(r)
    nm = len(mlist)
    try: 
        rho0 = np.load("mcdtens-U=0-B=0-m=%d-grid=%d-E=-1--0.5.npy" %(nm, nr))
    except: 
        print "!!!! cannot load rho0, continue without the actual data."
        rho0 = np.zeros((nr, nm, 4))
    print "rh0 loaded: max, min = ", np.max(rho0), np.min(rho0)
    rho = np.zeros((nr, nm, 4))    
    drho = np.zeros((nr, nm, 4))
    drhotot = np.zeros((nr, 4))
    Emax = -0.5
    Emin = -1.0
    emax = (np.abs(E-Emax)).argmin()
    emin = (np.abs(E-Emin)).argmin()
    h = E[emin + 1] - E[emin]
    for a in range (0,4):
        for m in range (0, nm):
            if (m == 0):
                for rr in [0, 1, 10, 20, 50]:
                   plot (E[emin:emax], mldostens[m, emin:emax, rr], 
                        label='r[%d]' % rr)
                   y_vals = mldostens[m, emin:emax, rr]
                   x_vals = E[emin:emax]
                   sum1 = sum (y_vals)
                   sum2 = sum (y_vals[1:-1])
                   sum3 = sum2 + 0.5 * y_vals[-1]
                   sum4 = sum2 + 0.5 * y_vals[0]
                   print "rr = ", rr, "sums: ", sum1, sum2, sum3, sum4
                legend()
                #show()
            for e1 in range (emin+1, emax-1):
                rho[:,m,a] += mldostens[m,e1,:] * (E[e1 + 1] - E[e1])
            if a == 1 or a == 3:
                rho[:,m,a] += mldostens[m,emin,:] * 0.5 * (E[emin+1]-E[emin])
            if a == 2 or a == 3:                
                rho[:,m,a] += mldostens[m,emax,:] * 0.5 * (E[emax-1]-E[emax-2])
   
    np.save("mcdtens-U=%g-B=%g-m=%d-grid=%d-E=%g-%g" %(info[0], info[1], nm, nr, Emin, Emax), rho)
    drho = rho - rho0
    theory = np.zeros((nr))
    for a in range(0,nr):
        theory[a] = (abs(Emin) - abs(Emax)) * (-0.05) #/ (2.0 * np.pi)
    np.save("mdrhos-U=%g-B=%g-m=%d-grid=%d-E=%g-%g" %(info[0], info[1], nm, nr, Emin, Emax), drho)
    if False:
        for m in range (0,nm):
            figure()
#            plot(r,rho0[:,m,0], label='rho0')
            plot(r,drho[:,m,0], label='drho no edges')
            plot(r,drho[:,m,1], label='drho with low edge')
            plot(r,drho[:,m,2], label='drho with high edge')
            plot(r,drho[:,m,3], label='drho both edges')
            title("Charge Density - -U=%g, m=%d" %(info[0], mlist[m]))
            legend()
    if True:
        for m in range (0,nm):
            for a in range (0,4):
                drhotot[:,a] += drho[:,m,a]
        figure()
        plot(r,theory, '--', label='Prediction')
        plot(r,drhotot[:,0], label='drho no edges')
        plot(r,drhotot[:,1], label='drho with low edge')
        plot(r,drhotot[:,2], label='drho with high edge')
        plot(r,drhotot[:,3], label='drho both edges')
        title("Total Induced Charge Density - -U=%g" %info[0])
#" %d m Channels" %(info[0], nm))
        legend()
        

    show()
    rho0tot = np.zeros((nr))
    figure()
    rho_expected = 1.0 / np.pi / 4.0 * (Emax * abs(Emax) - Emin * abs(Emin))
    plot([r[0], r[-1]], [rho_expected, rho_expected], 
         'k--', label="Expected value")
    for m in range (0,nm):
        rho0tot[:] += rho0[:,m,3]                 
    plot(r,rho0tot[:], label="total charge density, U = 0")
    legend()
    show()

    return 0

if __name__ == '__main__':
    mldostens = mldoscalc()
    mcdplotter(mldostens)
