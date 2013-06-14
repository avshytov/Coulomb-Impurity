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
    dostens = np.load("dostens-cn-U=%g-B=%g-ms=%d-N=%d-dm6.npy" 
                      %(info[0], info[1], nm, nr))
    cdtens = np.load("cdtens-cn-U=%g-B=%g-ms=%d-N=%d-dm6.npy" 
                     %(info[0], info[1],nm,nr))
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
    np.save("mldostens-cn-U=%g-ms=%d-B=%g-grid=%d-dm6" 
            %(info[0], nm, info[1], nr), mldostens)
    return mldostens

def mcdplotter(mldostens):
    r = np.load("rvec.npy")
    mlist = np.load("mlist.npy")
    E = np.load("Evec.npy")
    info = np.load("EMinfo.npy")
    nr = len(r)
    nm = len(mlist)
    Emax = -0.5
    Emin = -1.0
    Elims = [Emin, Emax]
    np.save("Elims",Elims)
    try: 
        rho0 = np.load("mcdtens-cn-U=0-B=0-m=%d-grid=%d-E=%g-%g-dm6.npy" 
                       %(nm, nr, Emin, Emax))
    except: 
        print "!!!! cannot load rho0, continue without the actual data."
        rho0 = np.zeros((nr, nm))
    print "rh0 loaded: max, min = ", np.max(rho0), np.min(rho0)
    rho = np.zeros((nr, nm))    
    drho = np.zeros((nr, nm))
    drhotot = np.zeros((nr))
    emax = (np.abs(E-Emax)).argmin()
    emin = (np.abs(E-Emin)).argmin()
    h = E[emin + 1] - E[emin]
    for m in range (0, nm):
        if False: #(m == 0):
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
        for e1 in range (emin+1, emax):
            rho[:,m] += mldostens[m,e1,:] * h
        rho[:,m] += mldostens[m,emin,:] * 0.5 * h
        rho[:,m] += mldostens[m,emax,:] * 0.5 * h

    np.save("mcdtens-cn-U=%g-B=%g-m=%d-grid=%d-E=%g-%g-dm6" 
            %(info[0], info[1], nm, nr, Emin, Emax), rho)
    drho = rho - rho0
    if abs(Emax - Emin) < 10e-6:
        theory = info[0]**2 / 2.0 / np.pi 
    else:
        theory = ((abs(Emin) - abs(Emax)) * -info[0]) / 2 / np.pi
    print "Theory =", theory
    np.save("mdrho-cn-U=%g-B=%g-m=%d-grid=%d-E=%g-%g-dm6" 
            %(info[0], info[1], nm, nr, Emin, Emax), drho)
    if False:
        for m in range (0,nm):
            figure()
            plot(r,drho[:,m], label='drho')
            title("Charge Density - U=%g, m=%d" %(-info[0], mlist[m]))
            legend()
    if False:
        for m in range (0,nm):
            drhotot[:] += drho[:,m]
        figure()
        plot([r[0],r[-1]], [theory, theory], '--', label='Prediction')
        plot(r,drhotot[:], label='drho')
        title("Total Induced Charge Density - U=%g" %-info[0])
#" %d m Channels" %(info[0], nm))
        legend()
        figure()
        plot(r,drhotot[:]/theory)
        title('Ratio Simulation/Theory, U=%g' %(-info[0]))
        show()

    if True:
        rho0tot = np.zeros((nr))
        figure()
        rho_expected = 1.0 / np.pi / 4.0 * (Emax * abs(Emax) - Emin * abs(Emin))
        plot([r[0], r[-1]], [rho_expected, rho_expected], 
             'k--', label="Expected value")
        for m in range (0,nm):
            rho0tot[:] += rho0[:,m]                 
        plot(r,rho0tot[:], label="total charge density, U = 0") ##
        legend()
        show()

    return 0

if __name__ == '__main__':
    mldostens = mldoscalc()
    mcdplotter(mldostens)
