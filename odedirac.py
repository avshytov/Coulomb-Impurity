import math
import numpy as np
import scipy 
from scipy import special



def _sgn(x): 
    if x > 0: return 1.0
    if x < 0: return -1.0
    return 0.0

sgn = np.vectorize(_sgn)

def odedos_m(E,r,U,m):

    def psireg(r_i,E): ### Not normalised here
        uu = special.jn(m, (E*r_i)) 
        ud = special.jn(m+1, (E*r_i))
        return uu,ud
    def psising(r_i, E): ### Not normalised here
        vu = special.yn(m, (np.abs(E)*r_i))
        vd = sgn(E) * special.yn(m+1, (np.abs(E)*r_i))
        return vu, vd
    
    chi_u = np.zeros((len(r), len(E)))
    chi_d = np.zeros(np.shape(chi_u))

    chi_u[0, :], chi_d[0, :] = psireg(r[0], E - U[0])

    nchi = np.sqrt(chi_u[0, :]**2 + chi_d[0, :]**2)
    chi_u[0, :] /= nchi
    chi_d[0, :] /= nchi

    if m >= 0:
        mu = m
    #    chi_u[0] = 1.0
    #    chi_d[0] = 0.0
    else:        
        mu = -m-1
    #    print 'mu =', mu
    #    chi_u[0] = 0.0
    #    chi_d[0] = 1.0
    def rhs(chi_u, chi_d, r, U): 
        f_u = (m-mu)/r * chi_u - (E - U)*chi_d
        f_d = (E - U)*chi_u - (1 + mu + m)/r * chi_d
        return f_u,f_d
    for i in range(1,len(r)):
        chi_up, chi_dp = chi_u[i-1, :], chi_d[i-1, :]
        h = r[i] - r[i-1]
        r1 = r[i-1]
        r2 = r[i-1] + 0.5 * h
        r3 = r2
        r4 = r[i]
        U1 = U[i-1]
        U2 = 0.5*(U[i-1]+U[i])
        U3 = U2
        U4 = U[i]
        k1u, k1d = rhs(chi_up, chi_dp, r1, U1)
        k2u, k2d = rhs(chi_up+k1u*0.5*h, chi_dp+0.5*k1d*h, r2, U2)
        k3u, k3d = rhs(chi_up+k2u*0.5*h, chi_dp+0.5*k2d*h, r3, U3)
        k4u, k4d = rhs(chi_up+k3u*h, chi_dp+k3d*h, r4 ,U4)
        
        chi_un = chi_up + h*(k1u + 2*k2u + 2*k3u + k4u)/6.0
        chi_dn = chi_dp + h*(k1d + 2*k2d + 2*k3d + k4d)/6.0

        chi_u[i, :], chi_d[i, :] = chi_un, chi_dn

    uu, ud = psireg(r[-1], E)
    vu, vd = psising(r[-1], E)
    D = uu * vd - vu * ud
    A = (chi_u[-1, :] * vd - chi_d[-1, :] * vu) / D
    B = (chi_d[-1, :] * uu - chi_u[-1, :] * ud) / D
    
    
    dos_m = ((r/r[-1])**(2*mu) * (np.abs(chi_u)**2 + np.abs(chi_d)**2).transpose()).transpose() 
    dos_m *= np.abs(E) / 4.0 / np.pi / (A**2 + B**2)
#rhomatch = A**2 * (abs(uu)**2+abs(ud)**2) + B**2 * (abs(vu)**2+abs(vd)**2)
    
    return dos_m

def doscalc(E,r,U,mlist):
    dos_tot = np.zeros((len(r), len(E)))
    for m in mlist:
        print m, E[0,], E[-1]
        dosm = odedos_m(E,r,U,m)
        dos_tot += dosm
    return 2.0 * dos_tot # We assume only m>=0 are included

def rhocalc(Emin, Emax, r, U, mlist):
    N_e = 1000
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
    Emin = -1.0
    Emax = -1e-4
    mlist = np.arange(0, 10.0, 1.0)
    r_0 = 1.0
    r = np.linspace(0.1, 50.0, 500)
    figure()
 
    for Z in [0.0, 0.1, 0.2, 0.3]:
        U = - Z / np.sqrt(r**2 + r_0**2)
        rho_0 = - (Emax**2 - Emin**2) / 4.0 / math.pi
        rho_rpa = -Z / 16.0 * r_0 / np.sqrt(r**2 + r_0**2)**3
        rho = rhocalc(Emin, Emax, r, U, mlist)
        plot (r, rho, label='Z = %g' % Z)
        plot (r, rho_0 + rho_rpa, label='Z = %g, rpa' % Z)
    legend()
    show()
    

if __name__ == '__main__':
    mlist = np.arange(0,10) 
    #plotode(1.0, mlist)
    test_rho()
   
