import numpy as np
import scipy
from scipy import linalg
import cmath
import math
from math import *
from pylab import *
from scipy import special
from ldos import *
import time

def makeDirac(r, U, m, B):
    N = len(r)
    #H = np.zeros((2*N,2*N), dtype=complex)
    P = np.zeros((2*N,2*N), dtype=complex)
    M = np.zeros((2*N,2*N), dtype=complex)
    V = np.zeros((2*N,2*N), dtype=complex)
    print "Calculating Momentum Channel:", m, "field B = ", B
    for y in range (0,N):
                if y == 0:
                    a = r[y+1] - r[y]
                    s = np.sqrt(r[1]/r[0])
                    P[0,1] = -1.0j / a * s
                    P[1,0] = 1.0j /a   * s
                else:
                    a = r[y] - r[y-1]
                    s = 1.0
                    if (y < N - 1):
                        s = np.sqrt(r[y + 1]/r[y])
                    s2 = 1.0 #np.sqrt(r[y-1]/r[y])
                    P[2*y+1,2*y]= 1.0j /a * s  ####FIX THIS
                    P[2*y,2*y+1]=  P[2*y+1,2*y].conj()
                    P[2*y,2*y-1]= 1.0j /a * s2 #* np.sqrt(r[y - 1]/r[y])
                    P[2*y-1,2*y]= P[2*y,2*y-1].conj() #-1.0j / a #* np.sqrt(r[y - 1]/r[y])
                # Calculating the m/r term. We split this 
                # into two steps, to allow for coupling to the 
                # left and right value of psi_1 --- AVS 
                if (y != N - 1): 
                    r_1 = r[y]
                    r_2 = r[y + 1]
                    y_next = y + 1
                    #y_next = y
                    #r_2 = r[y]
                else:
                    r_1 = r[-2]
                    r_2 = r[-1]
                    y_next = y
                r_mid = r_2 #(r_1 + 3.0*r_2)/4.0
                m_eff = m + 0.0 + B * r_mid**2 / 2.0
                M[2*y, 2*y_next + 1] += -0.0j * m_eff / r_mid
                #print "**", y, M[2 * y, 2*y_next + 1], 2*y, 2*y_next + 1
                M[2*y_next + 1, 2*y] = M[2*y, 2*y_next + 1].conj()
                #print M[2*y, 2*y_next + 1], M[2 * y_next + 1, 2 * y]
      
                if (y < N - 1): 
                    r_1 = r[y]
                    r_2 = r[y + 1]
                else:
                    r_1 = r[-2]
                    r_2 = r[-1]
                r_mid = r_1
                #r_mid = (3.0 * r_1 + r_2) / 4.0

                #if y < 6 :
                #    delta_m = np.sqrt(r_1) / (np.sqrt(r_1) + np.sqrt(r_2))
                #else:
                #    delta_m = 0.5
                #m_eff = m + delta_m + B * r_mid**2 / 2.0     
                m_eff = m + 0.0 + B * r_mid**2 / 2.0     
                #print m_eff, r_mid
                M[2*y, 2*y + 1] += -1.0j * m_eff / r_mid
                M[2*y + 1, 2*y] = M[2*y, 2*y + 1].conj()
                #print M[2 * y, 2*y + 1], 2 * y, 2 * y + 1
                #print M[2*y, 2*y + 1], M[2 * y + 1, 2 * y]
                V[2*y, 2*y] = U[y]
                V[2*y +1, 2*y +1] = U[y]
    
    H = P + M + V
    return H
    

def showWF(r, dr, E, vr, m):
    iplot = []
    for Eplot in [-0.24, -0.04]:
                Ei = list(enumerate (E))
                Ei.sort (lambda x, y: cmp(abs(x[1]- Eplot), abs(y[1]- Eplot)))
                iplot.append(Ei[0][0])
    for i in range (len(E)):
                u = vr[:,i]
                u_up = u[::2]
                u_down = u[1::2]
                if i in iplot:
                        print "Plotting state %d" %i
                        figure()
                        sqr = np.sqrt(r)
                        plot(r,u_up.real/sqr, label='psi up real')
                        plot(r,u_up.imag/sqr, label='psi up imag')
                        plot(r,u_down.real/sqr, label='psi down real')
                        plot(r,u_down.imag/sqr, label='psi down imag')
                        legend()
                        title('Momentum Channel %d, E = %g' % (m, E[i]))
    show()    
    
def findDensity(r, dr, E, vr):
    N = len(r)
    density = np.zeros((N, len(E))) 
    for i in range (len(E)):
                u = vr[:,i]
                u_up = u[::2]
                u_down = u[1::2]
                Dn = np.sum(abs(u_up)**2 + abs(u_down)**2)
                Dn -= 0.5 * (abs(u_up[0])**2 + abs(u_down[0])**2)
                Dn -= 0.5 * (abs(u_up[-1])**2 + abs(u_down[-1])**2)
                Cn = 1.0 / np.sqrt(Dn)
                density[:, i] =(abs(u_up)**2+abs(u_down)**2)/(2*np.pi*r*dr) * Cn
    return density

def solveDirac(r, U, m, B):
    N = len(r)
    dr = np.zeros((N))
    dr[1:] = r[1:] - r[0:-1]
    dr[0] = r[1] - r[0]
    #np.save("drvec", dr)
    
    H = makeDirac(r, U, m, B)
    print "diagonalising... "
    w, vr =  scipy.linalg.eigh(H)
    E = w
    if False:
                ea = -2
                eb = 2
                ivals = [t[0] for t in enumerate(w) if t[1] > ea and t[1] < eb]
                evals = [w[t] for t in ivals]
                if m == 0:
                    ens = list(evals)
                else:
                    ens.extend(evals)
    density = findDensity(r, dr, E, vr)
    #if m == 0: showWF (r, dr, E, vr, m)
    
    return E, density    

def dos_bg(mvals, Ev, r):
    bg = np.zeros((len(r), len(Ev)))
    
    for i in range(len(Ev)):
        s = np.zeros((len(r)))
        E = Ev[i]
        #print "E = ", Ev[i]
        for m in mvals:
            #print "add m = ", m
            jm = special.jn(m, E * r)
            jm1 = special.jn(m + 1, E * r)
            s += jm**2 + jm1**2
        s *= 1.0
        bg[:, i] = s / 2.0 / np.pi * abs(E)
    return bg 
    
def diracLDOS(Ev, r, U, mlist, B, gam):
    Nr = len(r)
    Nm = len(mlist)
    Ns = 2 * Nr
    
    N_e = len(Ev)
    LDOS = np.zeros ((Nr, N_e))
    A = 2.0/np.pi*gam**3 
    
    timestart1 = time.time() 
    for i_m, m in enumerate(mlist): 
        Bwvals = [(B, 1.0), (-B, 1.0)]
        if (abs(B) < 1e-2 * 1.0 / np.max(r)**2): # magnetic length too large
            Bwvals = [(0.0, 2.0)]
        
        for Bvalue, weight in Bwvals:
            E, d = solveDirac(r, U, m, Bvalue)
            for i in range(N_e):
                #X = A / (gam**2 + (Ev[i] - E)**2)**2       
                X = np.exp(-(Ev[i] - E)**2/gam**2) * np.sqrt(1.0/np.pi)/gam; 
                LDOS[:, i] += weight * np.dot(d, X)
    #LDOS += dos_bg(np.max(mlist) + 1, 100, Ev, r)
    #LDOS = dos_bg(0, 100, Ev, r)
    return LDOS

    
def find_rho(Ev, r, ldos, E_min, E_max):
    i_max = (np.abs(Ev - E_max)).argmin()
    i_min = (np.abs(Ev - E_min)).argmin()
    h = Ev[i_min + 1] - Ev[i_min]
    s = np.zeros(np.shape(r))
    s += np.sum(ldos[:, (i_min + 1):i_max], axis=1) 
    s += ldos[:, i_min] * 0.5
    s += ldos[:, i_max] * 0.5
    # It is very likely  that the window edge does not
    # coincide exactly with the grid point. In this case, 
    # one has to introduce the appropriate correction
    s +=  ldos[:, i_min] * (Ev[i_min] - E_min) / h
    s += -ldos[:, i_min] * (Ev[i_max] - E_max) / h
    print Ev[i_min], Ev[i_max]
    
    return s * h


def prepareLDOS(Ev, r, U0, U, mlist, B0, gam):
    fname = "ldos-U0=%g-N=%d-Mmax=%g-B=%g.npz" % (U0, len(r), np.max(mlist), B)
    try:
        data = np.load(fname)
        print "loaded data from", fname
        print "check Ev"
        #print 'Ev = ', Ev, 'data = ', data['Ev']
        #print 'diff = ', Ev - data['Ev']
        assert norm(Ev - data['Ev']) < 1e-8, "Ev are different"
        #print "check r"
        assert norm(r - data['r']) < 1e-8, "r vectors are different"
        #print "check U"
        assert norm(U - data['U']) < 1e-8, "U vectors are different"
        #print "check mlist"
        assert norm(mlist - data['mlist']) < 1e-8, "m lists are different"
        #print "check B"
        assert abs(B0 - data['B']) < 1e-4, "B values are different"
        print "check gamma"
        assert abs(gam - data['gam']) < 1e-6, "gamma values are different"
        print "OK, return precalculated ldos"
        return data['ldos']
    except:
        import traceback
        traceback.print_exc()
        print "cannot load data from file", fname, "recalculating"
        ldos = diracLDOS(Ev, r, U, mlist, B0, gam)
        np.savez(fname, ldos=ldos, Ev=Ev, r=r, U=U, mlist=mlist, B=B0, gam=gam)
        return ldos    
    
def graft(E_cut, Ev, r, ldos, r2, ldos2):    
    from scipy.interpolate import splev, splrep
    for i, E in enumerate(Ev):
        if abs(E) < E_cut:
           xvals = r2
           yvals = ldos2[:, i]
           spl = splrep (xvals, yvals)
           ldos[:, i] = splev(r, spl)
    
if __name__ == '__main__':
   ### dm terminates at element 6
   print "Running Coulomb potential."
   N = 1000
   rmin  = 0.01
   rmax  = 25.0
   Mmax  = 15 
   Mmax2 = 15 
   
   N2 = 500
   rmax2 = 200.0
   gam2 = np.pi * 0.8 / rmax2
   
   #gam = np.pi * 0.8 / rmax
   gam = np.pi * 0.8 / rmax
   
   #alpha = 0.2
   B0 = 0.0 #2.370
   r = zeros((N))
   U = zeros((N))
   
   E_cut = 0.1
   
   Ustr = 0.0
   r0 = 1.0
   
   #info = np.zeros((2))
   #info[0] = Ustr
   #info[1] = B0
   #np.save("EMinfo", info)
   
   mlist  = np.array(range(0, Mmax))
   mlist2 = np.array(range(0, Mmax2))
   
   r  = np.linspace(rmin, rmax,  N)
   r2 = np.linspace(rmin, rmax2, N2)
   U  = -Ustr /np.sqrt(r**2 + r0**2) #Coulomb
   U2 = -Ustr /np.sqrt(r2**2 + r0**2) #Coulomb
   #pot = -Ustr * np.exp(-alpha * r) #Exponential

   print "Momentum Channels:",  mlist
   Ev = np.linspace(-3.0, 3.0, 1001)
   ldos = prepareLDOS(Ev, r, Ustr, U, mlist, B0, gam)
   dos_0 = dos_bg(mlist, Ev, r)
   ldos2 = prepareLDOS(Ev, r2, Ustr, U2, mlist2, B0, gam2)
   
   graft(E_cut, Ev, r, ldos, r2, ldos2)
   
   figure()
   title ("DOS")
   for i_r, r_i in list(enumerate(r))[0:500:100]:
       plot(Ev, ldos[i_r, :], label='r = %g' % r_i)
       plot(Ev, dos_0[i_r, :], '--', label='BESSELS: r = %g' % r_i)
       ivals = [t for t in range(len(Ev)) if abs(Ev[t]) > 1e-2]
       xvals = np.array([np.log(abs(Ev[t])) for t in ivals])
       yvals = np.array([np.log(ldos[i_r, t]) for t in ivals])
       grad = np.polyfit(xvals, yvals, 1)
       print "r = ", r_i, "grad = ", grad
       alpha = grad[0]
       A  = grad[1]
       print "alpha = ", alpha, "A = ", A
       fitvals = alpha * xvals + A
       #alpha # + grad[1]
       print fitvals
       print np.shape(xvals), np.shape(fitvals)
       #plot(np.exp(xvals), np.exp(fitvals), '--', label='Fit: r = %g' % r_i)
   plot(Ev, np.abs(Ev) / 2.0 / np.pi, 'k--', label='linear')
   ivals = [t for t in range(len(Ev)) if Ev[t] > -1.0 and Ev[t] < 0.0]
   y_i = np.array([ldos[5, t] for t in ivals])
   t_i = np.array([np.abs(Ev[t]) / 2.0 / np.pi  for t in ivals])
   A_x = sum(y_i * t_i) / sum(t_i**2)
   print "A_x = ", A_x
   #plot (Ev, np.abs(Ev) / 2.0 / np.pi * A_x, 'r--', label='linear corrected')
   legend()

   
   
   #show()
   
   E_min = -1.0
   E_max = -0.1
   
   rho = find_rho(Ev, r, ldos, E_min, E_max)
   rho_B = find_rho(Ev, r, dos_0, E_min, E_max)
   
   #show()
   
   figure()
   title("no potential, density in the band [%g, %g]" % (E_min, E_max))
   plot (r, rho, label='simulation')
   plot (r, rho_B, label='BESSEL')
   rho_th = (abs(E_max) * E_max - E_min * abs(E_min)) / 4.0 / np.pi
   plot (r, 0*r + rho_th, 'k--', label='Theory')
   legend()
   #show()
   rho_0 = np.array(rho)
   def Uq_Coulomb(q):
       return 2.0 * np.pi / q * np.exp(-q*r0)
   from denscheck import polaris_generic
   #from rpakernel import kernel as rpa_kernel
   #K_RPA = rpa_kernel(r, mlist)
   rho_up   =  polaris_generic(r, E_max, Uq_Coulomb)
   rho_down =  polaris_generic(r, E_min, Uq_Coulomb)
   rho_RPA = rho_up - rho_down
   sgn = 1.0
   if (E_min < 0): sgn = -1
   rho_bg = -1.0 / 4.0 / np.pi * ((E_min - U) * np.abs(E_min - U) - E_min * abs(E_min)) 
   #rho_down + 1.0 / 4.0 / np.pi * U**2 * sgn
   
   results = []
  
   for U0 in [0.01, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]:
   #for U0 in [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4]:
       U  = -U0 /np.sqrt(r**2 + r0**2) #Coulomb
       U2 = -U0 /np.sqrt(r2**2 + r0**2) #Coulomb
       print "do: ", U0
       ldos  = prepareLDOS(Ev, r,  U0, U,  mlist,  B0, gam)
       ldos2 = prepareLDOS(Ev, r2, U0, U2, mlist2, B0, gam2)
       graft(E_cut, Ev, r, ldos, r2, ldos2)
       
       rho = find_rho(Ev, r, ldos, E_min, E_max)
       drho = rho - rho_0
       rho_bg = 1.0 / 4.0 / np.pi * ((E_min - U) * np.abs(E_min - U) - E_min * abs(E_min))
       drho_full = drho + rho_bg
       rho_TF = 1.0 / 4.0 / np.pi * ((E_max - U) * np.abs(E_max - U) - E_max * abs(E_max))
       results.append((U0, drho, drho_full, rho_TF))
       ivals = [t for t in range(len(r)) if r[t] > 1.0 and r[t] < 10.0]
       y_i = np.array([drho[t]/U0 for t in ivals])
       t_i = np.array([ rho_RPA[t] for t in ivals])
       A_x = sum(y_i * t_i) / sum(t_i**2)
       print "A_x = ", A_x
       if False:
          figure()
          title ("DOS: U0 = %g" % U0)
          for i_r, r_i in list(enumerate(r))[0:500:50]:
              plot(Ev, ldos[i_r, :], label='r = %g' % r_i)
          plot(Ev, np.abs(Ev) / 2.0 / np.pi, 'k--', label='linear')
          legend()
       
   if True:
      figure()
      title ("delta rho from the subband")
      for U0, drho, drho_full, rho_TF in results:
          plot (r, drho / U0, label='U0 = %g' % U0)
      plot (r, rho_RPA, 'k--', label='RPA response')
      legend()
   
   if True:
      figure()
      title ("delta rho from the band: log-log")
      for U0, drho, drho_full, rho_TF in results:
          loglog (r, abs(drho / U0), label='U0 = %g' % U0)
      loglog (r, abs(rho_RPA), 'k--', label='RPA response')
      legend()    
   #show() 
   
   if False:
      figure()
      title ("total delta rho")
      for U0, drho, drho_full, rho_TF in results:
          plot (r, drho_full, label="U_0 = %g" % U0)
          #plot (r, rho_TF, '--', label='TF, U_0 = %g' % U0)
      plot (r, rho_up, 'k--', label='RPA')
      legend()
   
   if True:
      figure()
      title ("total delta rho: log-log")
      for U0, drho, drho_full, rho_TF in results:
          loglog (r, np.abs(drho_full), label="U_0 = %g" % U0)
          #loglog (r, np.abs(rho_TF), '--', label="TF, U_0 = %g" % U0)
      loglog (r, np.abs(rho_up), 'k--', label='RPA')
      loglog (r, 1.0 / r, 'g--', label='1/r')
      loglog (r, 1.0/np.pi**2 / r**2, 'r--', label='1/r^2')
      legend()
   
   if False:
      figure()
      title ("nonlinear rho")
      for U0, drho, drho_full, rho_TF in results:
          plot (r, drho / U0 - rho_RPA, label='U0 = %g' % U0)
      legend()    

   
   if False:
      figure()
      title ("nonlinear rho: log-log")
      for U0, drho, drho_full, rho_TF in results:
          loglog (r, abs(drho / U0 - rho_RPA), label='U0 = %g' % U0)
      loglog (r, 1/r, 'k--', label='1/r')
      loglog (r, 1/r**2, 'r--', label='1/r^2')
      legend()    
      
   if False:
      for U0, drho, drho_full, rho_TF in results:
          figure()
          title("Band contribution, U0 = %g" % U0)
          plot (r, drho, label='sim')
          plot (r, rho_RPA * U0, label='RPA')
          plot (r, rho_TF - rho_bg, label='TF')
          legend()
   if True:
      for U0, drho, drho_full, rho_TF in results:
          figure()
          rho_bg = drho_full - drho
          title("Total density, U0 = %g" % U0)
          plot (r, drho_full, label='sim')
          plot (r, rho_up * U0, label='RPA + bg')
          #rho_bg = drho_full - drho
          plot (r, rho_TF, label='TF')
          legend()
   
   show()
   
   #E, dostens = DOS2(Ematp, Ematn, mlist ,r)
   #np.save("rvec",r)
   #np.save("mlist",mlist)
   #np.save("Evec",E)
   #np.save("potvec-cn-U=%g-grid=%g-r0=%g" %(Ustr, N, r0),pot)
