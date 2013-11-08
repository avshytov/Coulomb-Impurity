import math
import numpy as np
import scipy
import matplotlib
from matplotlib import pylab
from pylab import *
import fits

r = np.linspace(0.01,25,1000)
Us = np.load('Ustrengths.npy')
drhotots = np.load('drhomat.npy')
r0 = 1.0
rhoRPA = 1.0 / 16.0 / (r**2 + r0**2)**1.5
half = len(Us) / 2
modUs = np.array(Us[half::])
symrho=[]
asymrho=[]

for i in range (0,half):
    a = drhotots[(half-i-1),:] + drhotots[(half+i),:]
    a /= 2.0
    symrho.append(a)
    b = drhotots[(half+i),:] - drhotots[(half-i-1),:]
    b /= 2.0
    asymrho.append(b)
symrho=np.array(symrho)
asymrho=np.array(asymrho)
imin = 0
imax = 400
rmin = r[imin]
rmax = r[imax]
dr = r[imin+1]-r[imin]
Qsym =[]
Qasym =[]

if False:
    for i in range (0, len(Us)):
        figure()
        title('Total drho U0 = %g' %Us[i])
        plot(r,drhotots[i,:],label='sim')
        plot(r,Us[i]*rhoRPA,'k--',label='RPA')
        legend()

if False:
    for i in range (0,half):
        figure()
        title('symmetric drho U0=%g' %modUs[i])
        plot(r,symrho[i,:],label='sym')
        plot(r,drhotots[(half-i-1),:],'--',label='-ve U0')
        plot(r,drhotots[(half+i),:],'--',label='+ve U0')
        legend()

if False:
    for i in range (0,half):
        figure()
        title('asymmetric drho U0=%g' %modUs[i])
        plot(r,asymrho[i,:],label='asym')
        plot(r,drhotots[(half-i-1),:],'--',label='-ve U0')
        plot(r,drhotots[(half+i),:],'--',label='+ve U0')
        legend()

if False:
    figure()
    for i in range (0,half):
        title('Symmetric Scaled by U0^2')
        plot(r, symrho[i,:]/modUs[i]**2, label='U0=%g' %modUs[i])
        legend()

if False:
    figure()
    for i in range (0,half):
        title('Asymmetric Scaled by U0')
        plot(r, asymrho[i,:]/modUs[i], label='U0=%g' %modUs[i])
        legend()

if True:
    figure()
    for i in range (0,half):
        title('Symmetric Scaled by first point')
        plot(r, symrho[i,:]/symrho[i,0], label='U0=%g' %modUs[i])
        legend()

if True:
    figure()
    for i in range (0,half):
        title('Asymmetric Scaled by first point')
        plot(r, asymrho[i,:]/asymrho[i,0], label='U0=%g' %modUs[i])
        legend()

    figure()
    title('Asymmetric')
    for i in range (0,half):
        loglog(r, abs(asymrho[i,:]/asymrho[i,0]))
    legend()


if True:
    for k in range (0,half):
        Qs = 0.5 * dr * symrho[k,imin] * r[imin]
        Qs += 0.5 * dr * symrho[k,imax] * r[imax]
        Qas = 0.5 * dr * asymrho[k,imin] * r[imin]
        Qas += 0.5 * dr * asymrho[k,imax] * r[imax]
        for i in range (imin+1, imax):
            Qs += dr * symrho[k,i] * r[i]
            Qas += dr * asymrho[k,i] * r[i]
        Qs *= 2.0 * np.pi
        Qas *= 2.0 * np.pi
        Qsym.append(Qs)
        Qasym.append(Qas)
    Qtheory = (( 1.0 / np.sqrt(r0**2+rmin**2))
                - (1.0/ np.sqrt(r0**2 + rmax**2)))
    Qtheory *= modUs * np.pi / 8.0
    BS = np.pi/8.0 * modUs + 0.19 * modUs**3
    figure()
    title('Symmetric charge contribution')
    plot(modUs,Qsym, label='Sim')
    plot([modUs[0],modUs[-1]],[0,0],'--')
    symgrad = fits.linfit((Qsym/modUs**2), modUs**2)
    print 'Symmetric contribution =', symgrad[0], 'U0^4 + ', symgrad[1], 'U0^2'
    plot(modUs, (symgrad[0]*modUs**4+symgrad[1]*modUs**2), '--', label='fit')
    legend()

    figure()
    totQp = np.array(Qsym)+np.array(Qasym)
    totQn = np.array(Qasym)-np.array(Qsym)
    title('Antisymmetric charge contribution')
    plot(modUs,Qasym, label='Sim')
    plot(modUs,Qtheory,'--', label='Theory')
    plot(modUs,BS, 'k--', label='BS')
    cc, cl = fits.linfit(Qasym/modUs,modUs**2)
    print 'Asymmetric contribution =', cc, 'U0^3 + ', cl, 'U0'
    plot(modUs, cc*modUs**3 + cl*modUs, '--', label='fit')
    extra = Qasym - Qtheory
    legend()
    figure()

    title('additional asymmetric charge')
    loglog(modUs, abs(extra))
    loglog(modUs, 0.1*modUs**3, '--', label='U^3')
    loglog(modUs, 0.0001*modUs, '--', label='U')
    legend()
    figure()

    title('Total Charge with contirbutions')
    plot(modUs,Qsym,'--',label='Symmetric')
    plot(modUs,Qasym,'--',label='Antisymmetric')
    plot(modUs, totQp,label='Total Charge, positive U0')
    plot(modUs, totQn,label='Total Charge, negative U0')
    plot(modUs,Qtheory,'--',label='Theory')
    plot(modUs,BS, 'k--', label='BS')
    legend()
show()
