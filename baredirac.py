import math
from math import *
import cmath
import numpy
from numpy import *
import scipy
from scipy import special
from ldos import *
from diracsolver import *
from denscheck import polaris_generic

gam0 = 0.7

N=1000
rmin = 0.01
rmax = 25.0
Mmax = 15
gam = np.pi * gam0 / rmax
F = 1.01340014 + 0.05991426 / np.sqrt(N) + 7.12091516 / N #### Correction

N2=500
rmax2 = 200.0
Mmax2 = 15
gam2 = np.pi * gam0 / rmax2

B0=0.0

r = zeros((N))
U = zeros((N))

E_cut = 0.1

Ustr = 0.0
r0 = 1.0

mlist  = np.array(range(0, Mmax))
mlist2 = np.array(range(0, Mmax2))

r  = np.linspace(rmin, rmax,  N)
r2 = np.linspace(rmin, rmax2, N2)
U  = -Ustr /np.sqrt(r**2 + r0**2) #Coulomb   
U2 = -Ustr /np.sqrt(r2**2 + r0**2) #Coulomb

Ustrengths = [-0.1, 0.1]
def Uq_Coulomb(q):
    return 2.0 * np.pi / q * np.exp(-q*r0)

def Uvals(Beta, rvals):
    return -Beta /np.sqrt(rvals**2 + r0**2) #Coulomb

def highm(Ef, r, Mmax, U0):
    Jsum = np.zeros((len(r)))
    Efr = Ef + Uvals(U0,r)/2.0
    for i in range (-Mmax, Mmax+1):
        Jsum += 2.0 * (special.jn(i,(Ef*r)))**2
    drhohm = 2.0 - Jsum
    drhohm += (special.jn(-Mmax,(Ef*r)))**2
    drhohm -= (special.jn(1+Mmax,(Ef*r)))**2
    drhohm *= - abs(Ef) / 4.0 / np.pi * Uvals(U0,r)
    return drhohm

print "Momentum Channels:", mlist
print "U Strengths:", Ustrengths

Ev = np.linspace(-3.0, 3.0, 1001)
ldos = prepareLDOS(Ev, r, Ustr, U, mlist, B0, gam)
dos_0 = dos_bg(mlist, Ev, r)
ldos2 = prepareLDOS(Ev, r2, Ustr, U2, mlist2, B0, gam2)
####graft(E_cut, Ev, r, ldos, r2, ldos2)   

E_min = -1.0
E_max = -0.0
T = 3e-2

rho = find_rho(Ev, r, ldos, E_min, E_max)
rho_B = find_rho(Ev, r, dos_0, E_min, E_max)
rho_1 = getDensity(r, 0.0, 0.0*r, mlist, B0, E_min, E_max, T)
rho_0 = np.array(rho)
rho_wf_0 = np.array(rho_1)

rho_up   =  polaris_generic(r, E_max, Uq_Coulomb)
rho_down =  polaris_generic(r, E_min, Uq_Coulomb)
rho_RPA = rho_up - rho_down
sgn = 1.0
if (E_min < 0): sgn = -1
rho_bg = -1.0 / 4.0 / np.pi * ((E_min - U) * np.abs(E_min - U) - E_min * abs(E_min))

results = []

for U0 in Ustrengths:
    print "Calculating U0=%g" %U0
    U = Uvals(U0,r)
    U2 = Uvals(U0,r2)
    ldos  = prepareLDOS(Ev, r,  U0, U,  mlist,  B0, gam)
    ldos2 = prepareLDOS(Ev, r2, U0, U2, mlist2, B0, gam2)
    ####graft(E_cut, Ev, r, ldos, r2, ldos2)
    rho = find_rho(Ev, r, ldos, E_min, E_max)
    rho_1 = getDensity(r, U0, U, mlist, B0, E_min, E_max, T)
    rho_comp = 1.0 / (np.sqrt(r0**2 + r**2))**3 / 16.0
    drho = (rho - rho_0) * F
    rho_bg = 1.0 / 4.0 / np.pi * ((E_min-U) * np.abs(E_min-U) - E_min * abs(E_min))
    drho_full =  drho + rho_bg
    rho_TF = 1.0 / 4.0 / np.pi * ((E_max-U) * np.abs(E_max-U) - E_max * abs(E_max))
    drho_wf = rho_1 - rho_wf_0
    drhohm = highm(E_max,r,mlist[-1],U0) - highm(E_min,r,mlist[-1],U0)
    ###allrhos = [r,drho,drho_full,rho_TF,drho_wf,drhohm,rho_bg,U0*rho_up,U0*rho_down]
    ###np.save('allrhos-U0=%g-N=%d-Mmax=%d-B=%g-Emin=%g-Emax=%g' %(U0,N,np.max(mlist),B0,
    ###E_min,E_max), allrhos)
    results.append((U0, drho, drho_full, rho_TF, drho_wf))

