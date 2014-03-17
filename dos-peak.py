import sys
import numpy as np
from pylab import *
import odedirac
import os

#
#   A script to analyse DOS peak
#   python dos-peak.py file.npz Emin Emax
#   

if __name__ == '__main__':
    filename = sys.argv[1]
    Emin0 = float(sys.argv[2])
    Emax0 = float(sys.argv[3])
    data=np.load(filename)
    filename = filename[:-4]
    N0 = 1000
    rtest = 1.0
    threshold = 1e-9
    r = data['r']
    U = data['U']
    itest = abs(r-rtest).argmin()
    Mmax = 10
    mlist = np.arange(0.0, Mmax+0.5, 1.0)

    negEf = False
    fields = filename.split('-')
    for f in fields:
        if negEf == True:
            Ef = -float(f)
            negEf = False
        if f[:3] == 'Ef=':
            Ef = f[3:]
            if f[-1] == '=':
                negEf = True
        if f[:2] == 'Z=':
            Z = float(f[2:])
        if f[:3] == 'it=':
            it = int(f[3:])

    Ef = float(Ef)
    print 'Z =', Z, 'Ef =', Ef, 'it =', it
    testname = 'data/sim/denssim-Z=%g-N=500-Nfa=3.6-Ef=%g-it=%d.npz' %(Z,Ef,it+1)
    if os.path.isfile(testname):
        sys.exit('THIS IS NOT THE FINAL ITERATION!')
        
    Emax = Emax0
    Emin = Emin0

    Espace0 = np.linspace(Emin,Emax,N0)
    Espace = np.array(Espace0)
    dos0 = odedirac.doscalc(Espace0, r, U, mlist)
    dos = np.array(dos0)
    ipeak = dos0[itest,:].argmax()

    while (Emax - Emin)/N0 > threshold:
        if ipeak > 0:
            Emin = Espace[ipeak - 1]
        else: Emin = Espace[0]
        if ipeak < len(Espace) - 1:
            Emax = Espace[ipeak + 1]
        else: Emax = Espace[-1]
        Espace = np.linspace(Emin,Emax,N0)
        dos = odedirac.doscalc(Espace, r, U, mlist)
        ipeak = dos[itest,:].argmax()

    if abs(Emin-Emin0)<1e-7 or abs(Emax-Emax0)<1e-7:
        print 'CAREFUL, PEAK IS ON WINDOW EDGE!'

    Epeak = Espace[ipeak]
    height = dos[itest,ipeak]
    print "Peak position =", Epeak
    print "Peak height =", height
    chi_u, chi_d = odedirac.odepsi_m(Espace,r,U,0)

    def regions(E):
        return (E - U)**2 - (0.0 + 0.5)**2 / r**2

    Eplus = U + 0.5 / r
    Eminus = U - 0.5 / r

    figure()
    title('Wavefunctions for state E=%g' %Epeak)
    plot(r,chi_u[:,ipeak], label='chi up')
    plot(r, chi_d[:,ipeak], label='chi down')
    legend()

    figure()
    plot(r, U, label='Potential')
    plot(r,dos[:,ipeak], label='peak DOS(E=%g)' %Epeak)
    plot(r, 1e1*regions(Epeak), '--', label='SCALED Pr^2')
#    plot(r, Eplus, '--', label='Eplus')
#    plot(r, Eminus, '--', label='Eminus')
#    plot([r[0],r[-1]], [Epeak,Epeak], label='peak energy')
    axis([r[0], r[-1], 1.1*U[0], 1.1*dos[0,ipeak]])
    legend()

    filename = 'dospeaks-Z=%g.npz' %Z
    if os.path.isfile(filename):
        peakdata = np.load(filename)
        Efs = list(peakdata['Efs'])
        peaks = list(peakdata['peaks'])
        its = list(peakdata['its'])
        dosr = list(peakdata['dosr'])
        if Ef < Efs[0]:
            Efs.insert(0, Ef)
            peaks.insert(0, Epeak)
            its.insert(0, it)
            dosr.insert(0, dos[:,ipeak])
        elif Ef > Efs[-1]:
            Efs.append(Ef)
            peaks.append(Epeak)
            its.append(it)
            dosr.append(dos[:,ipeak])
        else:
            for i in range (0,len(Efs)):
                if abs(Efs[i]-Ef) < 1e-5:
                    peaks[i] = Epeak
                    its[i] = it
                    dosr[i] = dos[:,ipeak]
                elif (Efs[i]<Ef) and Efs[i+1]>Ef:
                    Efs.insert(i+1, Ef)
                    peaks.insert(i+1, Epeak)
                    its.insert(i+1, it)
                    dosr.insert(i+1, dos[:,ipeak])
        np.savez(filename[:-4], Efs=Efs, peaks=peaks, its=its, dosr=dosr)
    else:
        np.savez(filename[:-4], Efs=[Ef], peaks=[Epeak], its=[it], dosr=[dos[:,ipeak]])


    if False:
        figure()
        ra = 10.0
        rb = 40.0
        ivals = [t[0] for t in enumerate(r) if t[1] < rb and t[1] > ra]
        xvals = [r[t] for t in ivals]
        yvals = [(abs(U[t])) for t in ivals]
        grad = np.polyfit(log(xvals), log(yvals), 1)
        print "U tail exponent = ", grad[0]
        loglog(r, abs(U), label='U')
        loglog(xvals, (xvals**grad[0] * exp(grad[1])), label='fit')
        legend()

    if True:
        figure()
        plot(Espace0, dos0[itest, :], label='peak')
        plot(Espace, dos[itest,:], label='DOS at peak')
        legend()
    
    if False:
        from scipy import interpolate
        print 'Plotting 2D spatial peak LDOS'
        figure()
        X = np.linspace(-10.0, 10.0, 250)
        Y = np.array(X)
        spatdos = np.zeros((250,250))
        spl = interpolate.splrep(r,dos[:,ipeak])
        for n in range (len(Y)):
            d = np.sqrt(Y[n]**2 + X[:]**2)
            spatdos[:,n] = interpolate.splev(d[:], spl)
        pcolor(X,Y,spatdos,cmap='hot')
        colorbar()
        xlabel('x distance from impurity')
        ylabel('y distance from impurity')
    show()
