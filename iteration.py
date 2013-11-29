import numpy as np
from numpy import linalg
import math
from coulomb import coulombkernel
import scipy
from scipy import special
from scipy.interpolate import splev, splrep
from diracsolver import (makeDirac, solveDirac, dos_bg, diracLDOS, find_rho, getDensity,
                         prepareLDOS, bandDensity)
from pylab import *
import mkrpa2

def cheb_perm(r):
    """
        Construct the permutation which determines the 
        ordering of relaxation times tau.
        This is done recursively
    """
    if (r <= 0): return [ 0 ]
    if (r == 1):
        return [0, 1]
    prev = cheb_perm(r - 1)
    out = []
    nr = 2 * len(prev)
    for j in range(len(prev)):
        out.append(prev[j])
        out.append(nr - prev[j] - 1)
    return out;     

def cheb_tau(r, lam_max, lam_min):
    """ 
        Chebyshev set of relaxation parameters
        for a spectrum limited by eigenvalues lam_max and lam_min
    """
    i = 2**r; 
    tau_ch = np.zeros((i))
    av = (lam_max + lam_min)/2.0
    dif = (lam_max - lam_min) / 2.0
    for j in range(0, i):
        phi = (2. * j + 1.) * math.pi / 2.0 / float(i)
        tau_ch[j] = 1.0/(av + dif * math.cos(phi))
    return tau_ch; 

def make_tau_set(r, tau_max, tau_min):
    """
       Generate the set of relaxation parameters
    """
    #perm = [1, 8, 4, 5, 2, 7, 3, 6]
    #tau_max = 0.3
    #tau_min = 0.003; 
    taus = np.zeros((2**r))
    for i in range(len(taus)):
        taus[i] = tau_max * math.exp(math.log(tau_min/tau_max) * float(i)/(len(taus) - 1)) 
    #taus = [0.3, 0.2, 0.1, 0.05, 0.03, 0.01, 0.005, 0.003]
    #taus = cheb(3, lam_max, lam_min)
    perm = cheb_perm(r)
    print "perm: ", perm
    out = []
    for i in range(len(taus)):
        out.append(taus[perm[i] - 1])
    return out; 

N = 200
Mmax = 11
rexp = np.zeros((N))
rmin = 0.02
rmax = 40.0
for i in range (0,N):
    rexp[i] = rmin * math.exp(math.log(rmax/rmin)/(N - 1.0) * i)
    #r[i] = 0.1 * (i) + 0.01
r = np.linspace(rmin,rmax,N)
r_0 = 1.0
Ustr = 1.5
B0 = 0.0
Emin = -1.5
Emax = -0.0
Ev = np.linspace(-3.0, 3.0, 1001)
gam = 0.7 * np.pi / rmax
mlist = np.array(range(0, Mmax))
F = 1.01340014 + 0.05991426 / np.sqrt(N) + 7.12091516 / N #### Correction

#tau_u = 0.01
#tau_rho = 0.01
#0.3 0.2 0.1 0.05 0.03  0.01 0.005  0.003
#tau_u_set = [0.3, 0.003, 0.05, 0.03, 0.2, 0.005, 0.1, 0.01]
tau_u_set = make_tau_set(3, 0.2, 0.001)
#30.0, 1.0)
print "taus: ", tau_u_set
tau_rho_set = list(tau_u_set)
#tau_u_set = [0.3, 0.01, 0.1, 0.03]
#tau_rho_set = [0.3, 0.01, 0.1, 0.03]

U_0 = np.zeros((len(r)))
ldos_0 = prepareLDOS(Ev,r,0.0,U_0,mlist,0.0,gam)
rho_0 = bandDensity(r, U_0, mlist, 0.0, Emin, Emax, 0.01); 
if False:
   import pylab
   pylab.plot (r, rho_0 * 4.0 * 3.14159 / (Emax**2 - Emin**2))
   pylab.show()
#rho_0 = find_rho(Ev,r,ldos_0,Emin,Emax)

#np.savez('info-N=%d.npz' % (len(r)),Ustr=Ustr,B=B0,Emin=Emin,Emax=Emax,gam=gam,mlist=mlist,F=F,rho_0=rho_0)

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

def RPA_kernel_rho(U, r, kF):
    Q = mkrpa2.RPA_inter(r)  ## I cannot stand this nonsense! -- AVS
    if (kF * r.max() > 0.01):
        Q += mkrpa2.RPA_intra(r, kF)
    return np.dot(Q,U)

def highm(Ef,r,Mmax,U):
    Jsum = np.zeros((len(r)))
    Efr = Ef + U/2.0
    for i in range (-Mmax, Mmax+1):
        Jsum += special.jn(i, Ef*r)**2 + special.jn(i + 1, Ef*r)**2
        #Jsum += 2.0 * (special.jn(i,(Ef*r)))**2
    drhohm = 2.0 - Jsum
    #drhohm += (special.jn(-Mmax,(Ef*r)))**2
    #drhohm -= (special.jn(1+Mmax,(Ef*r)))**2
    drhohm *= - abs(Ef) / 4.0 / np.pi * U
    return drhohm

def seacontribution(r,rexp,rho,U,mlist,E_min,E_max):
    rho = highm(E_max,r,mlist[-1],U)
    rho -= highm(E_min,r,mlist[-1],U)
    rho -= U**2/4.0/np.pi
    rho += RPAresp(U,E_min,r,rexp)
    return rho

def RPAresp(U, Ef, r, rexp):
    Uexp = gridswap(r, rexp, U)
    # kF = abs(E_F) -- AVS
    # kernelrho was a silly name
    # Also, mind the spaces -- AVS
    rho_RPA_exp = RPA_kernel_rho(Uexp, rexp, abs(Ef))
    rho_RPA = gridswap(rexp, r, rho_RPA_exp)
    #Nend = 20
    #rho_RPA[-Nend:] = - abs(Ef) / 2.0 / math.pi * U[-Nend:]
    return rho_RPA

def rho_from_RPA(U, r, rexp, Ev):
    Nf = 4.0 
    alpha = 1.0
    Q = mkrpa2.RPA_inter(rexp)
    #for i in range(0, len(rexp)):
    #    Q[-i, :] = 0.0
    #    Q[-i, -i] = - 1.0 / 2.0 / math.pi * 0.1 
    #Q[-8:, :] = 0.0
    Uexp = gridswap(r, rexp, U)
    rho_RPA = np.dot(Q, Uexp)
    #rho_RPA[-20:] = 0.0
    rho = gridswap(rexp, r, rho_RPA)
    C = get_cou_mat(rexp)
    X = np.zeros((2 * len(rexp), 2*len(rexp)))
    X[0:len(rexp), len(rexp): ] = C * N * alpha
    X[len(rexp):,  0:len(rexp)] = Q
    
    #X = np.dot(C, Q) * 8.0 / math.pi 
    from scipy import linalg
    ev, w = linalg.eig(X, right=True)
    print ev - 1.0
    abs_ev = np.abs(ev - 1.0)
    abs_ev.sort()
    print abs_ev
    if False:
      import pylab 
      for i in [np.abs(ev - 1.0).argmax()]:
      #for i in range(len(ev)):
        #n0 = np.abs(ev).argmin()
        #if abs(ev[i] - 1.0) < 0.2:
            print i, ev[i], abs(ev[i])
            v_i = w[:, i]
            print "check: ", linalg.norm(np.dot(X, v_i) - ev[i] * v_i)
            rho0 = v_i[len(rexp):]
            u0   = v_i[:len(rexp)]
            print "rho: ", rho0
            print "u:", u0
            pylab.figure()
            pylab.plot(rexp, rho0.real, label='Re rho[%d]' % i)
            pylab.plot(rexp, rho0.imag, label='Im rho[%d]' % i)
            pylab.plot(rexp, u0.real, label='Re u[%d]' % i)
            pylab.plot(rexp, u0.imag, label='Im u[%d]' % i)
            pylab.plot(rexp, np.dot(Q, u0.real), label='Re Q*U0')
            pylab.plot(rexp, np.dot(Q, u0.imag), label='Im Q*U0')
            pylab.plot(rexp, np.dot(C, rho0.real), label='Re C*rho0')
            pylab.plot(rexp, np.dot(C, rho0.imag), label='Im C*rho0')
            pylab.title("ev[%d] = %g + i * %g abs = %g" % (i, ev[i].real, ev[i].imag, abs(ev[i])))
            pylab.legend()
      pylab.show()
    if False:
       import pylab
       pylab.figure()
       pylab.plot (rexp, rho_RPA, label='RPA')
       pylab.plot (r,    rho,     label='rho(r)')
       pylab.plot (rexp, -0.25/16.0 * r_0 / np.sqrt(rexp**2 + r_0**2)**3, label='rpa th')
       pylab.legend()
       pylab.figure()
       pylab.plot(rexp, Uexp, label='U')
       UU = -np.dot(C, rho_RPA)
       rho_u = 0.25 * r_0 / np.sqrt(rexp**2 + r_0**2)**3 / 2.0 / math.pi 
       pylab.plot(rexp, UU*8.0/3.14159, label='UU*...')
       pylab.plot(rexp, np.dot(C, rho_u), label='C*rho')
       pylab.legend()
       pylab.show()
    return rho * Nf * alpha

def rho_from_U(U,r,rexp,Ev):
    Nf = 4.0
    alpha = 1.0
    Rmax = 5.0
    #info = np.load('info-N=%d.npz' % (len(r)))
    #Ustr = info['Ustr']
    #B = info['B']
    #mlist = info['mlist']
    #E_min = info['Emin']
    #E_max = info['Emax']
    #F = info['F']
    #rho_0 = info['rho_0']
    #ldos = diracLDOS(Ev,r,U,mlist,B,gam)
    #rho = find_rho(Ev,r,ldos,E_min,E_max)
    rho = bandDensity(r, U, mlist, B0, Emin, Emax, 0.01)
    rho = F * (rho - rho_0)
    print "F = ", F
    if False:
        import pylab
        rho2 = rho + highm(Emax,r,mlist[-1],U)
        rho2 -= highm(Emin,r,mlist[-1],U)
        rho3 = rho2 - U**2/4.0/np.pi
        pylab.plot(r, rho, label='rho-Dirac, band')
        pylab.plot(r, rho2, label='band + high')
        pylab.plot(r, rho3, label='band + high + U2')
        #U1 = gridswap(r, rexp, U)
        #Q1 = mkrpa2.RPA_intra(rexp, abs(E_min))
        #if abs(E_max) * r.max() > 0.01:
        #   Q2 = mkrpa2.RPA_inter(rexp, abs(E_max))
        #else:
        #   Q2 = Q1 * 0.0; 
        #rho1 = np.dot(Q2 - Q1, U1)
        rho1 = RPAresp(U, Emax, r, rexp) - RPAresp(U, Emin, r, rexp)
        pylab.plot(r, rho1, label='Linear, band')
        pylab.legend()
        pylab.show()
    rho += seacontribution(r,rexp,rho,U,mlist,Emin,Emax) 
    imax = np.abs(r - Rmax).argmin()
    rho_rpa = RPAresp(U, Emax, r, rexp)
    rho[imax:] = rho_rpa[imax:]
    return rho * Nf * alpha

def gridswap(r1,r2,f1):
    spl = splrep(r1,f1)
    f2 = splev(r2,spl)
    return f2

def solve_coulomb(rho_U, U0, r, rexp, Ev, tau_u_set, tau_rho_set):
    info = np.load('info.npz')
    M = get_cou_mat(rexp)
    U = np.array(U0)
    rho = np.zeros((len(r)))
    it = 0
    j = 0
    Uerror = []
    rhoerror = []
    zero_U = True
    
    def mkFilename(it):
       fname = 'rhoandpot-U0=%g-B=%g-N=%g-Mmax=%g-it=%d.npz' %(
                   Ustr, B0, len(r), np.max(mlist), it)
       return fname
    
    if it > 0:
       data = np.load( mkFilename( it ) )
       rho = data['rho']
       U = data['U']
    while True:
        tau_U   = tau_u_set  [it % len(tau_u_set)]
        tau_rho = tau_rho_set[it % len(tau_rho_set)]
        it += 1
        rho = gridswap(r,rexp,rho)
        U1 = np.dot(M, rho)
        U1 = gridswap(rexp,r,U1)
        rho = gridswap(rexp,r,rho)
        U1 += U0
        if zero_U: U1 -= U1[-1];
        err_U = linalg.norm(U - U1)
        U += tau_U * (U1 - U)
        if zero_U: U -= U[-1];
        rho1 = rho_U(U, r, rexp, Ev)
        err_rho = linalg.norm(rho1 - rho)
        rho += tau_rho * (rho1 - rho)
        print "U error", err_U, "rho error", err_rho, "it", it
        np.savez(mkFilename( it ), rho=rho, U=U, rho1=rho1, U1=U1, r=r)
        Uerror.append(err_U)
        rhoerror.append(err_rho)
        if it % 10 == 0:
            np.savez('Errors-U0=%g-B=%g-N=%d.npz' 
                     %(info['Ustr'],info['B'],len(r)),
                     erU=Uerror, errho=rhoerror)
        if False:
            if it % 5 == 0:
                j+=1
                if j < 8:
                    figure(0)
                    title('Potential')
                    plot(r,U, label='it=%d' %it)
                    legend()
                    figure(1)
                    title('Rho')
                    plot(r,rho, label='it=%d' %it)
                    legend()
                else:
                    figure(0)
                    title('Potential')
                    plot(r,U,'--', label='it=%d' %it)
                    legend()
                    figure(1)
                    title('Rho')
                    plot(r,rho,'--', label='it=%d' %it)
                    legend()
                if j == 13:
                    j = 0
                    figure(0)
                    plot(r,U0, label='U0')
                    legend()
                    figure(2)
                    plot(Uerror, label='U error')
                    plot(rhoerror, label='rho error')
                    legend()
                    show()
        if (err_U < (1e-7)) and (err_rho < (1e-7)):
            break

    return U, rho

Nf = 4.0 ###
alpha = 2.5

def rho_minus12(U,r,a):
    return -U / r / 2.0 / np.pi**2 * Nf * alpha

U0 = Ustr * (r**2 + r_0**2)**(-0.5)

U, rho = solve_coulomb(rho_from_U, U0, r, rexp, Ev, tau_u_set, tau_rho_set) 
#U, rho = solve_coulomb(rho_from_RPA, U0, r, rexp, Ev, tau_u_set, tau_rho_set) 

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
