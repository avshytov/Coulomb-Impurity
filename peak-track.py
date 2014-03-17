import numpy as np
import matplotlib
from pylab import *

rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman']})

Econv = 1e9 * 1.05457173e-34 * 1e6 / 1.60217657e-19 ###eV
if Econv != 1.0:
    print 'Energy conversion:', Econv

Zs = [2, 2.5, 3, 3.5, 4, 4.5, 5]
cols = ['r', 'g', 'b','y','m','c','k']

botend = [True, True, True, True, False, False, False, False]
topend = [True, True, True, True, False, False, False, False]
### Last Ef where there is a classically
### forbidden region
botfor = [1.0, -0.05, -0.13, -0.22, -0.3, -100.0, -100.0]
topfor = [1.0, -0.02, -0.02, 0.01, 0.02, 0.04, 0.05]

for i in range (len(Zs)):
    Z = Zs[i]
    data = np.load('dospeaks-Z=%g.npz' %Z)
    Efs = data['Efs'] * Econv
    peaks = data['peaks'] * Econv
    ibot = abs(np.array(Efs) - botfor[i]*Econv).argmin()
    itop = abs(np.array(Efs) - topfor[i]*Econv).argmin()

    boundEfs = Efs[ibot:itop]
    boundpeaks = peaks[ibot:itop]

    plot(Efs, peaks, cols[i], ls='dotted')
    plot(boundEfs, boundpeaks, cols[i], label='Z=%g' %Z)
    plot(Efs, peaks, cols[i]+' .')
    if botend[i]:
        plot(Efs[0], peaks[0], cols[i]+' o')
    if topend[i]:
        plot(Efs[-1], peaks[-1], cols[i]+' o')
    
    text(Efs[0]-0.012, peaks[0]+0.002, '$Z$=%g' %Z , fontsize=14)

plot([-0.4, 0.4],[0.0, 0.0], 'k--')
plot([-1, 1], [-1, 1], 'k--')
axis([-0.26, 0.08, -0.035, 0.102])
xlabel('E$_{F}$ (eV)')
ylabel('E$_{0}$ (eV)')
title('DOS Peaks vs Fermi Energy')

#savefig('plots/Ef-Peak-Graph.pdf', paperstyle='a3')

show()
