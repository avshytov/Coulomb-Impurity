import numpy as np
import math
import mkrpa2
from scipy import integrate, interpolate, special, linalg

                                                             
def do_corr_intra(r, kF):
    """
       Calculate the intraband kernel corrections resulting from
       r > r[-1]
    """
    Npow = 3
    Q1 = np.zeros((len(r), Npow + 1))
    #if True: 
    #   pylab.figure()
    integrate_all = False
    dr0 = min(0.2/kF, 0.2)
    Rmax = 2.0 * max(r) + 100.0 / kF;
    Gfun = mkrpa2.mk_intra_spline(kF, Rmax)
    Rstar = r[-1]
    ximin = Rstar / Rmax; 
    print "ximin = ", ximin
    
    for i in range (0, len(r)):
        ri = r[i]
        #
        # Integrate from r[-1] to infinity, 
        #         
        # r = Rstar / xi, 0 < xi < 1
        #
        # Assume U ~ r^{-p}, i.e., U(r) = (Rstar/r)**p = xi^p
        #
        for p in range(Npow + 1):
            def Fp(xi):
                return Gfun(r[i], Rstar/xi) / xi**(3 - p) * Rstar**2
            Ip, epsp = integrate.quad(Fp, ximin, 1.0, limit=10000)
            print Ip, epsp
            if False:
               xivals = np.arange(ximin, 1.0, 1e-4)
               fvals = np.vectorize(Fp)(xivals)
               import pylab as pl
               pl.loglog(xivals, np.abs(fvals))
               pl.title('r = %g p = %g' % (ri, p))
               pl.show()
            Q1[i, p] = Ip
        print "intra-corr: ", i, Ip, epsp    
    return - Q1 * 2.0 * math.pi 

def RPA_corr_intra(r, kF, label=''):
    
    #integrate_all = False
    """
       Attempt to load the intraband kernel from file, 
       calculate if unavailable
    """
    Rmin = r.min()
    Rmax = r.max()
    N = len(r)
    fname = "data/rpakernel-intra-corr-kF=%g-Rmin=%g-Rmax=%g-N=%d-%s.dat.npz" % (kF, Rmin, Rmax, N, label)
    try: 
        data = np.load(fname)
        print "Intraband kernel correction loaded from", fname
        assert linalg.norm(r - data['r']) < 1e-6, "r grids are different"
        return data['Q']
    except:
        import traceback
        traceback.print_exc()
        print "cannot load correction data from", fname, "; recalculating"
        Q = do_corr_intra(r, kF)
        np.savez(fname, r=r, Q=Q, kF=kF)
        return Q
           
def do_corr_inter(r):
    """
       Calculate the interband kernel corrections 
       from r > r[-1]
    """
    integrate_all = True
    Qs = mkrpa2.mk_inter_spline()
    N = len(r)
    Npow = 3
    Q = np.zeros ((N, Npow + 1))
    Rstar = r[-1]
    for i in range (0, len(r)):
        ri = r[i]
        #
        # Integrate from r[-1] to infinity, 
        #         
        # r = Rstar / xi, 0 < xi < 1
        #
        # Assume U ~ r^{-p}, i.e., U(r) = (Rstar/r)**p = xi^p
        #
        for p in range(Npow + 1):
            def Fp(xi):
                return Qs(ri/Rstar * xi, 1.0) * xi**p
            Ip, epsp = integrate.quad(Fp, 0.0, 1.0, limit=1000)
            Q[i, p] = Ip
        print "inter-corr:", i
    return Q / Rstar

def RPA_corr_inter(r, label=''):
    """
       Attempt to load interband kernel from file, 
       recalculate if the data is not available
    """
    Rmin = r.min()
    Rmax = r.max()
    N = len(r)
    fname = "data/rpakernel-inter-corr-Rmin=%g-Rmax=%g-N=%g-%s.dat.npz" % (Rmin, Rmax, N, label)
    try: 
        data = np.load(fname)
        print "Interband kernel correction loaded from", fname
        assert linalg.norm(r - data['r']) < 1e-6, "r grids are different"
        return data['Q']
    except:
        import traceback
        traceback.print_exc()
        print "cannot load correction data from", fname, "; recalculating"
        Q = do_corr_inter(r)
        np.savez(fname, r=r, Q=Q)
        return Q

if __name__ == '__main__':
    
   rmin = 0.01
   rmax = 50.0
   N = 500
   r = rmin * np.exp(math.log(rmax/rmin)/(N - 1.0) * np.arange(0, N, 1.0))
   #r = np.arange(0.01, 10.0, 0.01)
   kF = 0.3
   Q_inter = mkrpa2.RPA_inter(r, 'exp')
   Q_intra  = mkrpa2.RPA_intra(r, kF, 'exp')
   Q_inter_corr = RPA_corr_inter(r, 'exp')
   Q_intra_corr = RPA_corr_intra(r, kF, 'exp')
   print "at 0: inter corr:", Q_inter_corr[0, :] * 2.0 * math.pi
   print "at 0: intra corr:", Q_intra_corr[0, :] * 2.0 * math.pi 
   print "sum(Q_intra) = ", sum(Q_intra[0, :]) * 2.0 * math.pi 
   r_0 = 1.0
   U = 1.0 / (r**2 + r_0**2)**0.5
   def Uq(q):
       return 2.0 * math.pi / (q + 1e-8) * math.exp(-q*r_0)
   
   def rho_direct(x):
       s = 1.0; #max(1.0, kF * x)
       def f(q):
           return Uq(q/s) * mkrpa2.Pi_full(q/s, kF) * special.jn(0, q * x / s) * q / s**2
       I1, eps1 = integrate.quad(f, 0, 2.0, limit=500)
       I2, eps2 = integrate.quad(f, 2.0, 4.0, limit=500)
       I3, eps3 = integrate.quad(f, 4.0, np.inf, limit=500)
       I = I1 + I2 + I3
       return - I / 2.0 / math.pi
   
   rho_RPA = np.dot(Q_inter, U)
   
   
   C_RPA = 0.0
   for i in range(len(r) - 1):
       y = (rho_RPA[i]*r[i] + rho_RPA[i + 1]*r[i + 1]) / 2.0
       C_RPA += y * (r[i + 1] - r[i])
   C_RPA += rho_RPA[0] * r[0]**2/2.0
   C_RPA += rho_RPA[-1] * r[-1]**2
   C_RPA *= 2.0 * math.pi
   print "Charge: ", C_RPA, "should be", -math.pi/8.0
   
   corr_0  = Q_inter_corr[:, 1] * U[-1] 
   corr_kf = Q_intra_corr[:, 1] * U[-1] 
   
   C_corr = 0.0
   for i in range(len(r) - 1):
       y = (corr_0[i]*r[i] + corr_0[i + 1]*r[i + 1]) / 2.0
       C_corr += y * (r[i + 1] - r[i])
   C_corr += corr_0[0] * r[0]**2/2.0
   C_corr += corr_0[-1] * r[-1]**2
   C_corr *= 2.0 * math.pi
   print "Charge correction: ", C_corr, ", total", C_RPA + C_corr, ", should be", -math.pi/8.0
   
   
   rho_TH = -1.0 / 16.0 * r_0 / (r**2 + r_0**2)**1.5
   rho_intra = np.dot(Q_intra, U)
   rho_th_intra = -U * abs(kF) / 2.0 / math.pi
   rho_d = np.vectorize(rho_direct)(r)
   import pylab
   pylab.plot (r, rho_RPA, label='Numerics-kf=0')
   pylab.plot (r, rho_RPA + corr_0, label='Numerics-kf=0 + corr')
   pylab.plot (r, rho_RPA + rho_intra, label='Numerics-kf=%g' % kF)
   pylab.plot (r, rho_RPA + rho_intra + corr_0 + corr_kf, 
               label='Numerics-kf=%g + corr' % kF)
   pylab.plot (r, corr_0, label='corr kf = 0')
   pylab.plot (r, corr_kf, label='corr_intra')
   pylab.plot (r, rho_TH,  label='Exact-kf=0')
   pylab.plot(r, -U*abs(kF)/2.0/math.pi, label='TF kf=%g' % kF)
   pylab.plot(r, rho_d, label='Direct integration')
   pylab.legend()
   #pylab.ylim(-1.0, 1.0)
   pylab.figure()
   pylab.loglog (r, np.abs(rho_RPA), label='Numerics')
   pylab.loglog (r, np.abs(rho_RPA + corr_0), label='Numerics + corr')
   pylab.loglog (r, np.abs(corr_0), label='corr_kF=0')
   pylab.loglog (r, np.abs(rho_RPA + rho_intra), label='Numerics-kf=%g' % kF)
   pylab.loglog (r, np.abs(rho_RPA + rho_intra + corr_0 + corr_kf), 
                    label='Numerics-kf=%g + corr' % kF)
   pylab.loglog (r, np.abs(rho_TH),  label='Exact')
   pylab.plot(r, abs(U)*abs(kF)/2.0/math.pi, label='TF kf=%g' % kF)
   pylab.loglog(r, np.abs(rho_d), label='Direct integration')
   pylab.legend()
   pylab.show()
   
