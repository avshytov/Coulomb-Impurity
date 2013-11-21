import numpy as np
from numpy import linalg
import math
from coulomb import coulombkernel
import scipy
from scipy import special
from scipy.interpolate import splev, splrep
from diracsolver import (makeDirac, solveDirac, dos_bg, diracLDOS, find_rho, getDensity,
                         prepareLDOS)

N = 1000
Mmax = 11
rexp = np.zeros((N))
rmin = 0.01
rmax = 25.0
for i in range (0,N):
    rexp[i] = rmin * math.exp(math.log(rmax/rmin)/(N - 1.0) * i)
    #r[i] = 0.1 * (i) + 0.01
r = np.linspace(rmin,rmax,N)
r_0 = 1.0
Ustr = 0.25
B0 = 0.0
Emin = -1.0
Emax = -0.0
Ev = np.linspace(-3.0, 3.0, 1001)
gam = 0.7 * np.pi / rmax
mlist = np.array(range(0, Mmax))
F = 1.01340014 + 0.05991426 / np.sqrt(N) + 7.12091516 / N #### Correction

tau_u = 0.1
tau_rho = 0.1

U_0 = np.zeros((len(r)))
ldos_0 = prepareLDOS(Ev,r,0.0,U_0,mlist,0.0,gam)
rho_0 = find_rho(Ev,r,ldos_0,Emin,Emax)

np.savez('info.npz',Ustr=Ustr,B=B0,Emin=Emin,Emax=Emax,gam=gam,mlist=mlist,F=F,rho_0=rho_0)

def get_cou_mat(r):
    N=len(r)
    fname = "coumat-N=%g.npz" %N
    try:
        coumat = np.load(fname)
        print "Coulomb Kernel loaded"
        print "Check r"
        assert linalg.norm(r - coumat['r'])<1e-8, "r vectors match"
        return coumat['kernel']
    except:
        import traceback
        traceback.print_exc()
        print "cannot load", fname, ": recalculating"
        coumat = coulombkernel(r)
        np.savez(fname, kernel=coumat, r=r)
        return coumat

def highm(Ef,r,Mmax,U):
    Jsum = np.zeros((len(r)))
    Efr = Ef + U/2.0
    for i in range (-Mmax, Mmax+1):
        Jsum += 2.0 * (special.jn(i,(Ef*r)))**2
    drhohm = 2.0 - Jsum
    drhohm += (special.jn(-Mmax,(Ef*r)))**2
    drhohm -= (special.jn(1+Mmax,(Ef*r)))**2
    drhohm *= - abs(Ef) / 4.0 / np.pi * U
    return drhohm

def seacontribution(r,rho,U,mlist,E_min,E_max):
    rho += highm(E_max,r,mlist[-1],U)
    rho -= highm(E_min,r,mlist[-1],U)
    rho -= U**2/4.0/np.pi
    rho += RPAresp(U,E_min)
    return rho

def RPAresp(U,Ef):
    return -abs(Ef)/2.0/np.pi*U

def rho_from_U(U,r,Ev):
    info = np.load('info.npz')
    Ustr = info['Ustr']
    B = info['B']
    mlist = info['mlist']
    E_min = info['Emin']
    E_max = info['Emax']
    F = info['F']
    rho_0 = info['rho_0']
    ldos = diracLDOS(Ev,r,U,mlist,B,gam)
    rho = find_rho(Ev,r,ldos,E_min,E_max)
    rho = F * (rho - rho_0) + rho_0
    return seacontribution(r,rho,U,mlist,E_min,E_max)

def gridswap(r1,r2,f1):
    spl = splrep(r1,f1)
    f2 = splev(r2,spl)
    return f2

def solve_coulomb(rho_U, U0, r, rexp, Ev, tau_U, tau_rho):
    info = np.load('info.npz')
    M = get_cou_mat(rexp)
    U = np.array(U0)
    rho = np.zeros((len(r)))
    it = 0
    j = 0
    Uerror = []
    rhoerror = []
    while True:
        it+=1
        fname = 'rhoandpot-U0=%g-B=%g-N=%g-Mmax=%g-tauU=%g-taurho=%g-it=%d' %(
            info['Ustr'],info['B'],len(r),np.max(info['mlist']),tau_U,tau_rho,it)
        rho = gridswap(r,rexp,rho)
        U1 = np.dot(M, rho)
        U1 = gridswap(rexp,r,U1)
        rho = gridswap(rexp,r,rho)
        U1 += U0
        err_U = linalg.norm(U - U1)
        U += tau_U * (U1- U)
        rho1 = rho_U(U, r, Ev)
        err_rho = linalg.norm(rho1 - rho)
        rho += tau_rho * (rho1 - rho)
        print "U error", err_U, "rho error", err_rho, "it", it
        np.savez(fname, rho=rho, U=U, rho1=rho1, U1=U1)
        Uerror.append(err_U)
        rhoerror.append(err_rho)
        if it % 10 == 0:
            np.savez('Errors-U0=%g-B=%g-N=%d-tauU=%g-taurho=%g.npz' 
                     %(info['Ustr'],info['B'],len(r),tau_U,rau_rho),
                     erU=Uerror, errho=rhoerror)
#        if True:
#            if it % 1 == 0:
#                j+=1
#                if j < 8:
#                    figure(0)
#                    title('Potential')
#                    plot(r,U, label='it=%d' %it)
#                    legend()
#                    figure(1)
#                    title('Rho')
#                    plot(r,rho, label='it=%d' %it)
#                    legend()
#                else:
#                    figure(0)
#                    title('Potential')
#                    plot(r,U,'--', label='it=%d' %it)
#                    legend()
#                    figure(1)
#                    title('Rho')
#                    plot(r,rho,'--', label='it=%d' %it)
#                    legend()
#                if j == 13:
#                    j = 0
#                    figure(0)
#                    plot(r,U0, label='U0')
#                    legend()
#                    figure(2)
#                    plot(Uerror, label='U error')
#                    plot(rhoerror, label='rho error')
#                    legend()
#                    show()
        if (err_U < (1e-7)) and (err_rho < (1e-7)):
            break

    return U, rho

Nf = 4.0 ###
alpha = 2.5

def rho_minus12(U,r,a):
    return -U / r / 2.0 / np.pi**2 * Nf * alpha

U0 = Ustr * (r**2 + r_0**2)**(-0.5)

U, rho = solve_coulomb(rho_from_U, U0, r, rexp, Ev, tau_u, tau_rho) 

ra = 5.0
rb = 15.0

ivals = [t[0] for t in enumerate(r) if t[1] < rb and t[1] > ra]
xvals = [r[t] for t in ivals]
yvals = [(abs(U[t])) for t in ivals]

grad = np.polyfit(log(xvals), log(yvals), 1)
print grad
fitvals = exp(grad[1] + grad[0]*log(xvals))

#figure()
#plot (r, U)
#figure()
#loglog (r, abs(U), label='U')
#loglog (r, 1.0/r, label='1/r')
#loglog(r, 1.0/r/r, label='1/r^2')
#loglog(xvals, abs(fitvals), label='Fit: alpha = %g' % (-grad[0]))
#legend()
#show()
