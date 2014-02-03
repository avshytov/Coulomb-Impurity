import math
import numpy as np
import scipy 
from scipy import special
from scipy import interpolate
import util 
from util import * 

def _sgn(x): 
    if x > 0: return 1.0
    if x < 0: return -1.0
    return 0.0

sgn = np.vectorize(_sgn)

def psireg(r_i, Ex, m): ### Not normalised here                                          
    uu = special.jn(m, (np.abs(Ex)*r_i))
    ud = special.jn(m+1, (np.abs(Ex)*r_i)) * sgn(Ex)
    return uu, ud

def psising(r_i, Ex, m): ### Not normalised here                                         
    vu = special.yn(m, (np.abs(Ex)*r_i))
    vd = sgn(Ex) * special.yn(m+1, (np.abs(Ex)*r_i))
    return vu, vd


#
# We write the wave function as
#
#  psi_u =   (r/R)^mu chi_u
#  psi_d = i (r/R)^mu chi_d
#
#  where mu = m for m >= 0, and -m-1 otherwise. R == r[-1] here
#  This way, we can start with finite solution at r=0 and integrate
#  it to r = R. After that, we match it to wave function of free particle:
#  
#   psi = A * psi_reg(E, r, m) + B * psi_sing(E, r, m)
#
#  We then normalise the wave function to unit intensity: A^2 + B^2 = 1. 
#

def odepsi_m(E, r, U, m):
    
    chi_u = np.zeros((len(r), len(E)))
    chi_d = np.zeros(np.shape(chi_u))

    chi_u[0, :], chi_d[0, :] = psireg(r[0], E - U[0], m)

    nchi = np.abs(chi_u[0, :]) + np.abs(chi_d[0, :]) # this does not involve
                                                     # squares of small numbers
    np.sqrt(chi_u[0, :]**2 + chi_d[0, :]**2)         # unlike prev version

    chi_u[0, :] /= nchi
    chi_d[0, :] /= nchi

    if m >= 0:
        Ku = 0
        Kd = -1.0 - 2.0 * m
        mu = m
    else:        
        mu = -m-1
        Ku = 1.0 + 2.0 * m
        Kd = 0
        
    #
    #  Right-hand side of the ODE: chi'(r) = rhs(r)
    #
        
    def rhs(chi_u, chi_d, r, U):
        f_u = Ku / r * chi_u - (E - U)*chi_d
        f_d = (E - U)*chi_u  + Kd / r * chi_d
        #f_u = (m-mu)/r * chi_u - (E - U)*chi_d
        #f_d = (E - U)*chi_u - (1 + mu + m)/r * chi_d
        return f_u,f_d
    Us = interpolate.splrep(r, U)
    
    #
    # RK4 fourth order method
    #
    def rk4step(chi_u, chi_d, r_p, r_n, h):
        r1 = r_p
        r2 = r_p + 0.5 * h
        r3 = r2
        r4 = r_n
     
        #rs = np.array([r1, r2, r3, r4])
        Ui = interpolate.splev(np.array([r1, r2, r3, r4]), Us)

        U1 = Ui[0]
        U2 = Ui[1]
        U3 = Ui[2]
        U4 = Ui[3]

        k1u, k1d = rhs(chi_u, chi_d, r1, U1)
        k2u, k2d = rhs(chi_u+k1u*0.5*h, chi_d+0.5*k1d*h, r2, U2)
        k3u, k3d = rhs(chi_u+k2u*0.5*h, chi_d+0.5*k2d*h, r3, U3)
        k4u, k4d = rhs(chi_u+k3u*h, chi_d+k3d*h, r4 ,U4)

        chi_un = chi_u + h*(k1u + 2*k2u + 2*k3u + k4u)/6.0
        chi_dn = chi_d + h*(k1d + 2*k2d + 2*k3d + k4d)/6.0
        return chi_un, chi_dn

    for i in range(1,len(r)):
        chi_un, chi_dn = chi_u[i-1, :], chi_d[i-1, :]
        h = r[i] - r[i-1]
        dxi = max(abs(E)) * h
        # dxi is the phase change over this grid step
        # We split the step into smaller sub-steps if dxi 
        # is too large
        dxi0 = 0.1
        n_steps = int(dxi / dxi0) + 1
        dr = h / n_steps
        # Go over sub-steps
        for i_step in range(n_steps):
            r_p = r[i-1] + i_step * dr
            r_n = r[i-1] + (i_step + 1) * dr
            chi_un, chi_dn = rk4step(chi_un, chi_dn, r_p, r_n, dr)

        chi_u[i, :], chi_d[i, :] = chi_un, chi_dn
    uu, ud = psireg(r[-1], E, m)
    vu, vd = psising(r[-1], E, m)
    D = uu * vd - vu * ud
    A = (chi_u[-1, :] * vd - chi_d[-1, :] * vu) / D
    B = (chi_d[-1, :] * uu - chi_u[-1, :] * ud) / D
    r_mu =  (r/r[-1])**(mu)
    #
    # The goal of the crazy expression below is to multiply
    # chi by the corresponding values of r^mu. In numpy, multiplying
    # a matrix chi(r, E) by vector r would not work: we have to 
    # transpose the matrix first. 
    #
    psi_u =  (r_mu * chi_u.transpose()).transpose()
    psi_d =  (r_mu * chi_d.transpose()).transpose()
    psi_u /= np.sqrt(A**2 + B**2)
    psi_d /= np.sqrt(A**2 + B**2)

    return psi_u, psi_d

def odedos_m(E, r, U, m):
    psi_u, psi_d = odepsi_m(E, r, U, m)
    dos_m = np.abs(psi_u)**2 + np.abs(psi_d)**2
    dos_m *= np.abs(E) / 4.0 / np.pi
    return dos_m; 


def doscalc(E,r,U,mlist):
    dos_tot = np.zeros((len(r), len(E)))
    for m in mlist:
        print m, E[0,], E[-1]
        dosm = odedos_m(E, r, U, m)
        dos_tot += dosm
    return 2.0 * dos_tot # We assume only m>=0 are included

def rhocalc(Emin, Emax, r, U, mlist):
    N_e = 2000
    Es = np.linspace(Emin, Emax, N_e)
    dos = doscalc(Es, r, U, mlist)
    weight = np.zeros(len(Es))
    weight[:] = 1.0
    weight[0] = weight[-1] = 0.5
    dE = Es[1] - Es[0]
    rho = np.dot(dos, weight) * dE 
    return rho

def plotode(Emax, mlist):
    r = np.arange(0.01, 40.0, 0.05)
    U = 0.0 * r
    E = 1.0
    m = 0
    Es = np.arange(-Emax, Emax, 0.01)
    #rho = odedos_m(E,r,U,m)
    dos = doscalc(Es,r,U,mlist)
    for r_i in [0.1, 1.0, 2.0, 5.0, 10.0, 20.0]:
        i = np.abs(r - r_i).argmin()
        plot(Es, dos[i, :], label='r = %g' % r[i])
    plot (Es, np.abs(Es) / 2.0/math.pi, '--k', label='Expected')
    legend()
    title ('DOS')
    
    #rho_0 = special.jn(m, E*r - r*U)**2 + special.jn(m+1, E*r - r*U)**2
    #rho_0 *= abs(E) / 4.0 / np.pi
    #title('Charge Density')
    #plot(r, rho, label='rho from ode m=%d' %m)
    #plot(r, rhomatch, label='rhomatch m=%d' %m)
    #plot(r, rho_0, 'r--', label='bessel')
    #plot(r,rho/rhomatch)
    #legend()
    #figure()
    #title('density of states')
    #plot(r, dos, label='dos for E=&g')
    
    show() 


def test_rho():
#    import pylab as pl
    Emin = -1.0
    Emax = -1e-4
    mlist = np.arange(0, 10.0, 1.0)
    r_0 = 1.0
    r = np.linspace(0.1, 50.0, 500)
#    pl.figure()
 
    for Z in [0.0, 0.1, 0.2, 0.3]:
        U = - Z / np.sqrt(r**2 + r_0**2)
        rho_0 = - (Emax**2 - Emin**2) / 4.0 / math.pi
        rho_rpa = -Z / 16.0 * r_0 / np.sqrt(r**2 + r_0**2)**3
        rho = rhocalc(Emin, Emax, r, U, mlist)
#        pl.plot (r, rho, label='Z = %g' % Z)
#        pl.plot (r, rho_0 + rho_rpa, label='Z = %g, rpa' % Z)
#    pl.legend()
#    pl.show()
    

if __name__ == '__main__':
    mlist = np.arange(0,10) 
    #plotode(1.0, mlist)
    test_rho()
   
