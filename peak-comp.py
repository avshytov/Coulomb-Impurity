import numpy as np
from pylab import *
import odedirac
import lorentzfit as lfit

N0 = 1000
### Value of Z, lower and upper limits of EF for analysis
Zlist = [(3, -0.13, -0.02),
         (3.5, -0.18, 0.01),
         (4, -0.3, 0.02),
         (4.5, -0.3, 0.04),
         (5, -0.35, 0.05)]

window = 0.3
windows = [0.01, 0.005, 0.005, 0.005]
mlist = np.arange(0,0.5,1.0)
rtest = 1.0

def sgn(x):
    if x < 0: return -1.0
    else: return 1.0

for l in range (0, len(Zlist)):
    Z = Zlist[l][0]
    Elow = Zlist[l][1]
    Ehigh = Zlist[l][2]
    data = np.load("dospeaks-Z=%g.npz" %Z)
    ilow = abs(Elow-data['Efs']).argmin()
    ihigh = abs(Ehigh-data['Efs']).argmin() + 1

    allEfs = data['Efs']
    allpeaks = data['peaks']
    allits = data['its']
    allheights = data['dosr']
    Efs = allEfs[ilow:ihigh]
    peaks = allpeaks[ilow:ihigh]
    its = allits[ilow:ihigh]
    heights = allheights[ilow:ihigh]

    ### Finding Width ###
    widths = []
    dosr = []
    hts = []
    hfs = []

    for i in range (0, len(its)):
        Epeak = peaks[i]
        Ef = Efs[i]
        it = its[i]
        simdata = np.load('data/sim/denssim-Z=%g-N=500-Nfa=3.6-Ef=%g-it=%d.npz'
                              %(Z, Ef, it))
        U = simdata['U']
        r = simdata['r']
        itest = abs(r - rtest).argmin()
        height = heights[i][itest]

        normE = np.linspace(0.1, 4.0, 1000)
        trueE = normE * Epeak
        normdos = odedirac.doscalc(trueE, r, U, mlist)[itest,:] / height

        print 'Z =', Z, 'Ef =', Ef
        
        bnds = []
        bnds.append((0.9, 1.1)) # peak amplitude
        bnds.append((0.9, 1.1)) # peak position
        bnds.append((0.1, 1.5)) # peak width
        bnds.append((0.2, 1.2)) # envelope E scaling
        bnds.append((0.01, 2.0*normdos[-1])) # bg height
        bnds.append((0.5, 5.0)) # bg E scaling
        bnds = tuple(bnds)

        start = [[1.0, 1.0, 1.0, 0.5, 1.5*normdos[-1], 1.5]]

        fit, params = lfit.peakfit(normE, normdos, 1, bnds, start)
        hf = params[0][0]
        pf = params[0][1]
        wf = params[0][2]
        eo = params[0][3]
        bgh = params[0][4]
        bgo = params[0][5]
        widths.append(wf)
        hfs.append(hf)
        print 'tabulated Epeak', Epeak, 'Fitted Epeak scaling factor', pf
        print 'Lorentzian height', hf, 'width', wf
        print 'bg height', bgh, 'bg Escaling', bgo
        print 'envelope scaling', eo

        bg = bgh*np.tanh(normE/bgo)
        purelor = hf / ((normE-pf)**2 + wf**2) * wf / np.pi
        lorenv = purelor * np.tanh(normE/eo)
        peakonly = normdos - bg

        if False:
            figure()
            title('Ef = %g' %Ef)
            xlabel('E/Epeak')
            ylabel('DOS/DOSpeak')
            xlim([-1.0, 5.0])
            ylim([0.0, 1.3])
            plot(normE, normdos, label='scaled DOS')
            plot(normE, fit, 'k', label='fit')
            plot(normE, bg, '--', label='bg')
            plot(normE, purelor, label='Lorentzian')
            plot(normE, np.tanh(normE/eo), '--', label='envelope')
            plot(normE, peakonly, '--', label='DOS-bg')

        legend()
        if i < 7:
            figure(2*l)
#            plot(normE, normdos, label='dos Ef=%g' %Ef)
            plot(normE, peakonly/peakonly.max(), label='dos-bg Ef=%g' %Ef)
            figure(2*l+1)
            plot(normE, lorenv/lorenv.max(), label='Lorentzian Ef=%g' %Ef)
        elif i >= 7 and i < 14:
            figure(2*l)
#            plot(normE, normdos,'--',  label='dos Ef=%g' %Ef)
            plot(normE, peakonly/peakonly.max(), '--', label='dos-bg Ef=%g' %Ef)
            figure(2*l+1)
            plot(normE, lorenv/lorenv.max(), '--', label='Lorentzian Ef=%g' %Ef)
        else:
            figure(2*l)
#            plot(normE, normdos, '*', label='dos Ef=%g' %Ef)
            plot(normE, peakonly/peakonly.max(), '*', label='dos-bg Ef=%g' %Ef)
            figure(2*l+1)
            plot(normE, lorenv/lorenv.max(), '*', label='Lorentzian Ef=%g' %Ef)

        figure(2*l)
        title('Normalised DOS - Background, Z=%g' %Z)
        legend()
        figure(2*l+1)
        title('Normalised fitted Lorentzian x envelope, Z=%g' %Z)
        legend()
        dosr.append(heights[i])
        hts.append(heights[i][itest])

    widths = np.array(widths)
    hts = np.array(hts)
    hfs = np.array(hfs)

    figure(20)
    title('Areas, Epeak^0.5*Height*widths')
    plot(Efs, np.sqrt(abs(peaks))*hts*widths, label='Z=%g' %Z)
    plot(Efs, np.sqrt(abs(peaks))*hts*widths, '.')
    legend()

    figure()
    title('normalised DOS at peak Z=%g' %Z)
    for t in range (0, len(its)):
        if t < 7:
            semilogy(r, (dosr[t]*abs(peaks[t])**0.5), 
                     label='Ef = %g peak = %g' %(Efs[t], peaks[t]))
        elif t >= 7 and t < 14:
            semilogy(r, (dosr[t]*abs(peaks[t])**0.5), '--', 
                     label='Ef = %g peak = %g' %(Efs[t], peaks[t]))
        else:
            semilogy(r, (dosr[t]*abs(peaks[t])**0.5), '*', 
                     label='Ef = %g peak=%g' %(Efs[t], peaks[t]))
    legend()

show()
