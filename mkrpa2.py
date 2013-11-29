import numpy as np
import math
from scipy import integrate, interpolate, special, linalg


def Pi_intra (q):
    """
       Regular part of the polarization operator. The full 
       formula for polarization operator is
           Pi = q/16 + Pi_intra * kF / (2.0 * math.pi)
       At large values of q, this quantity decays which results
       in regular behaviour at r = r'
    """
    if (q <= 2): return (1.0 - math.pi/8.0 * q) 
    sq = math.sqrt(1.0 - 4.0/q/q)
    asin = math.asin(2.0/q)
    return (1.0 - 0.5 * sq - q/4.0*asin) 

def Pi_full(q, kF):
    return Pi_intra(q/kF) * kF / 2.0 / math.pi + q/16.0

def Pi_intra2 (q):
    """
       Pi_intra with asymptotic behaviour, 2/3q^2 subtracted. 
       I have changed 1/q^2 to 1/(q^2 + 1), because (a) it is finite 
       at q=0, and (b)the integra for the latter can be easily computed
    """
    return Pi_intra(q) - 2.0 / 3.0 / (q**2 + 1.0)


def show_pi():
   qvals = np.arange(0.001, 5.0, 0.01)
   qvals1 = np.arange(1.00, 5.0, 0.01)
   
   Pi_intra_vals = np.vectorize(Pi_intra)(qvals)
   Pi_full_vals  = np.vectorize(lambda q: Pi_full(q, 1.0))(qvals) 
   q2inv = np.vectorize(lambda q: 2.0/3.0/(q**2 + 1.0))(qvals)
   import pylab
   pylab.plot (qvals, Pi_intra_vals, label='Pi_intra')
   pylab.plot (qvals, Pi_full_vals, label='Pi_full')
   pylab.plot (qvals, q2inv, label='2/3q**2')
   pylab.plot(qvals, Pi_intra_vals - q2inv, label='diff')
   pylab.legend()
   pylab.show()


def F_intra2(r1, r2, kF):
    #print r1, r2, kF
    """
        Integration routine --- for internal use in Q_intra
    """
    s = np.max([r1 * kF, r2 * kF, 1.0])
    #
    # First, handle the difference between Pi_intra 
    # and 2/3 1/(q^2 + 1)
    #
    def f(x):
        J1 = special.jn(0, x * kF * r1 / s)
        J2 = special.jn(0, x * kF * r2 / s)
        return Pi_intra2(kF * x / s) * J1 * J2 * x / s**2
    #print r1, r2
    I, eps = integrate.quad(f, 0.0, np.Inf, limit=300)
    
    #
    # Now add the missing contribution:
    #   The integral for 1/(q^2 + 1) can be found as 
    #   Int q dq /(q^2 + 1) J_0(qr)  ~ K_0(r)
    #   With two Bessel functions, we can apply addition theorem, 
    #   which gives K_0(|r - r'|), averaged over angles
    #
    def f2(theta):
        r = math.sqrt(r1**2 + r2**2 + 2 * r1 * r2 * math.cos(theta))
        return special.kn(0, kF * r)
    I2, eps2 = integrate.quad(f2, 0.0, math.pi)
    I2 *= 1.0 / math.pi * 2.0 / 3.0 
    
    #print "In F2: ", I, I2
    return (I + I2) * kF**3 / (4.0 * math.pi**2)

def G2_intra(r, kF): # xi == kF * r
    print "r, KF = ", r, kF
    #print r1, r2, kF
    """
        Integration routine --- for internal use in Q_intra
    """
    s = np.max([kF * r, 1.0])
    #
    # First, handle the difference between Pi_intra 
    # and 2/3 1/(q^2 + 1)
    #
    def f(x):
        J1 = special.jn(0, x * kF * r / s)
        return Pi_intra2(x * kF / s) * J1 * x / s**2
    #print r1, r2
    I, eps = integrate.quad(f, 0.0, np.Inf, limit=300)
    
    #
    # Now add the missing contribution:
    #   The integral for 1/(q^2 + 1) can be found as 
    #   Int q dq /(q^2 + 1) J_0(qr)  ~ K_0(r)
    #   With two Bessel functions, we can apply addition theorem, 
    #   which gives K_0(|r - r'|), averaged over angles
    #
    #def f2(theta):
    #    r = math.sqrt(r1**2 + r2**2 + 2 * r1 * r2 * math.cos(theta))
    #    return special.kn(0, kF * r)
    #I2, eps2 = integrate.quad(f2, 0.0, math.pi)
    I2 = special.kn(0, kF * r)
    I2 *=  2.0 / 3.0 
    
    #print "In F2: ", I, I2
    return (I + I2) * kF**3 / (4.0 * math.pi**2)

def F_intra(r1, r2, kF):
    """
        Integration routine --- for internal use in Q_intra
        Here, we calculate intraband contribution in more direct
        way. This routine is slow and less reliable than F_intra2
        Provided mostly for testing F_intra2
    """
    s = max(kF * r1, kF * r2, 1.0)
    def f(x):
        J1 = special.jn(0, x * kF * r1 / s)
        J2 = special.jn(0, x * kF * r2 / s)
        return Pi_intra(kF * x / s) * J1 * J2 * x / s**2
    I, eps = integrate.quad(f, 0.0, np.Inf, limit=10000)
    I2 = 0.0
    print "In F1:", I
    return (I + I2) * kF**3 / (4.0 * math.pi**2)

    
    
def mk_intra_spline(kF, rmax):
    xvals = np.arange(0.001, max(kF * rmax, 1.0), 0.03)
    yvals = np.vectorize(lambda r: G2_intra(r, kF))(xvals)
    print xvals, yvals
    spl = interpolate.splrep(xvals, yvals)
    def F_intra(r1, r2):
        def f_theta(theta):
            R = math.sqrt(r1**2 + r2**2 + 2.0 * r1 * r2 * math.cos(theta))
            return interpolate.splev(R, spl, der=0)
        I, eps = integrate.quad(f_theta, 0, math.pi) 
        return I  / math.pi
    return F_intra

def testG(r1vals, r2vals, kF):
    import pylab
    pylab.figure()
    Gfun = mk_intra_spline(kF, max(r1vals) + max(r2vals))
    for r1 in r1vals:
        pylab.figure()
        F2vals = np.vectorize(lambda r2: F_intra2(r1, r2, kF))(r2vals)
        Gvals  = np.vectorize(lambda r2: Gfun(r1, r2))(r2vals)
        pylab.plot (r2vals, F2vals, label='F2')
        pylab.plot (r2vals, Gvals, label='G')
        pylab.title('r1 = %g' % r1)
    pylab.show()

def testF(r1, r2, kF):
    """
        Test F_intra2 by comparing it to F_intra
        The results should match to 1e-6 - 1e-8
    """
    F1 = F_intra(r1, r2, kF)
    F2 = F_intra2(r1, r2, kF)
    Gfun = mk_intra_spline(kF, 10.0)
    G = Gfun(r1, r2)
    print "r1, 2 = ", r1, r2, "kF = ", kF
    print "F1 = ", F1, "F2 = ", F2, "diff = ", F1 - F2, "rel: ", (F1 - F2)/(F1 + F2)*0.5
    print "G = ", G, "diff = ", G - F1, "rel:", (G - F1)/(G + F1) * 0.5

def do_RPA_intra(r, kF):
    """
       Calculate the intraband kernel
    """
    Q1 = np.zeros((len(r), len(r)))
    dr = np.zeros((len(r)))
    dr[1:] = np.diff(r)
    dr[0] = r[1] - r[0]
    #if True: 
    #   pylab.figure()
    integrate_all = False
    dr0 = 0.5
    Gfun = mk_intra_spline(kF, max(r) * 2.0)
    for i in range(len(r)):
        print "Qintra:", i
        for j in range(0, len(r) - 1):
            #if i == j: continue
           # if (j  == len(r) - 1):
           #     r1 = r[-2]
           #     r2 = r[-1]
           # else:
            r1 = r[j]
            r2 = r[j + 1]
            rc = (r1 + r2) / 2.0
            dr = r2 - r1
            #print r1, r2
            
            nj = int (dr / dr0) + 1
            #print r[j], nj, dr, dr0, dr/dr0
            Fk = np.zeros((nj,))
            xk = np.linspace(r1, r2, nj + 2)[1:-1]
            #Fk = np.zeros ((nj + 3,))
            #xk = np.linspace(r1, r2, nj + 3)
            for k in range(nj):
                Fk[k] = Gfun(r[i], xk[k]) * xk[k]
                #Fk[k] = F_intra2(r[i], xk[k], kF) * xk[k]  
            dxk = dr / nj
            
            I1 = sum(Fk) * dxk
            I2 = sum(Fk * (xk - rc)) * dxk
            
            Q1[i, j]     += 0.5 * I1 - I2 / dr
            Q1[i, j + 1] += 0.5 * I1 + I2 / dr            
            if False: #integrate_all:
               def G1(x):
                   #print "G1:", r[i], x, kF
                   return F_intra2(r[i], x, kF) * x
               def G2(x):
                   #print "G2:", r[i], x, kF
                   return F_intra2(r[i], x, kF) * x * (x - rc)
               I1, eps1 = integrate.quad(G1, r1, r2)
               I2, eps2 = integrate.quad(G2, r1, r2)
               #if (j == len(r) - 1):
               #    
               #else:
               Q1[i, j]     += 0.5 * I1 - I2 / dr
               Q1[i, j + 1] += 0.5 * I1 + I2 / dr
            #else: if False:
            #   Fij = F_intra2(r[i], rc, kF)
            #   Q1[i, j]     += 0.5 * Fij * r[j] * dr
            #   Q1[i, j + 1] += 0.5 * Fij * r[j] * dr
            #Fij = F_intra2(r[i], r[j], kF)
            #Q1[i, j] = Fij * r[j] * dr[j]
            #Q1[j, i] = Fij * r[i] * dr[i]
        Q1[i, 0] += F_intra2(r[i], 0.0, kF)* r[0]**2 / 2.0
        
        print "sum: ", sum(Q1[i, :]) * 4.0 * math.pi**2 # must be one
        if False and i % 20 == 0:
           import pylab
           pylab.figure()
           pylab.plot(r, Q1[i, :], label='r = %g' % r[i])
           pylab.show()
    return - Q1 * 2.0 * math.pi

def RPA_intra(r, kF):
    
    integrate_all = False
    """
       Attempt to load the intraband kernel from file, 
       calculate if unavailable
    """
    Rmin = r.min()
    Rmax = r.max()
    N = len(r)
    fname = "rpakernel-intra-kF=%g-Rmin=%g-Rmax=%g-N=%d.dat.npz" % (kF, Rmin, Rmax, N)
    try: 
        data = np.load(fname)
        print "Intraband kernel loaded from", fname
        assert linalg.norm(r - data['r']) < 1e-6, "r grids are different"
        return data['Q']
    except:
        import traceback
        traceback.print_exc()
        print "cannot load data from", fname, "; recalculating"
        Q = do_RPA_intra(r, kF)
        np.savez(fname, r=r, Q=Q, kF=kF)
        return Q
        
        

def F_inter(r1, r2, eps):
    """
       Calculate the interband kernel
    """
    def f(theta):
        R = math.sqrt(eps**2 + r1*r1 + r2*r2 - 2 * r1 * r2 * math.cos(theta))
        return (1.0 - 3.0*eps**2/R**2) / R**3
    I, epsI = integrate.quad(f, 0.0, math.pi)
    return I / 16.0  / math.pi

def mk_inter_spline():
    """
       Construct spline interpolation of F_inter
       as a function of r1/r2 (r1 < r2), 
       and a function that reuses this interpolation
    """
    def F(x):
        return F_inter(x, 1.0, 1e-2)
    xvals = np.arange(0.0001, 1.0001, 0.001)
    yvals = np.vectorize(F)(xvals)
    spl = interpolate.splrep(xvals, yvals)
    def Q_spline(r1, r2):
        if (r1 <= r2):
            x = r1 / r2
            r = r2
        else:
           x = r2 / r1
           r = r1
        y = interpolate.splev(x, spl, der=0)
        return y / r**3
    return Q_spline
   
def do_RPA_inter(r):
    """
       Calculate the interband kernel
    """
    integrate_all = False
    Qs = mk_inter_spline()
    N = len(r)
    Q = np.zeros ((N, N))
    for i in range (0, N):
        print "Qinter: ", i
        ri = r[i]
        # 
        # Integrate from r[0] to ri. Introduce 
        # dimensionless variable v = r / ri, 0 < v < 1
        #
        for j in range (0, i):
            v1 = r[j] / ri
            v2 = r[j + 1] / ri
            dv = v2 - v1
            vc = (v1 + v2)/2.0
            def f1(v):
                return Qs(ri, ri * v)
            def f2(v):
                return Qs(ri, ri * v) * (v - vc) 
            if integrate_all or (abs(i - j) < 10) or (i < 20):
                 I1, eps1 = integrate.quad(f1, v1, v2)
                 I2, eps2 = integrate.quad(f2, v1, v2)
            else:
                 I1  = Qs(ri, ri * vc) * dv
                 qs1 = Qs(ri, ri * v1)
                 qs2 = Qs(ri, ri * v2)
                 I2  = (qs2 - qs1) * dv**2 / 12.0
            #print "<", i, j, I1 * ri, I2/dv * ri
            Q[i, j]     += ( I1 / 2.0 - I2 / dv ) * ri**2 * v1 
            Q[i, j + 1] += ( I1 / 2.0 + I2 / dv ) * ri**2 * v2
        
        #
        # Now integrate from r' to r[-1]: r = ri / u,   0 < u < 1
        #
        for j in range (i + 1, N):
            u1 = ri/r[j]
            u2 = ri/r[j - 1]
            du = u2 - u1
            r1 = r[j - 1]
            r2 = r[j]
            uc = (u1 + u2)/2.0
            def f3(u):
                return Qs(ri, ri/u) / u**3
            def f4(u):
                return Qs(ri, ri/u) / u**3 * (u - uc) 
            if integrate_all or (abs(i - j) < 10) or (i < 20):
                I3, eps3 = integrate.quad(f3, u1, u2)
                I4, eps4 = integrate.quad(f4, u1, u2)
            else:
                I3  = Qs(ri, ri/uc) / uc**3 * du
                qs1 = Qs(ri, ri/u1) / u1**3  
                qs2 = Qs(ri, ri/u2) / u2**3 
                I4 = (qs2 - qs1)*du**2 / 12.0  
            #print ">", i, j, I3 * ri**2, I4/du * ri**2
            Q[i, j - 1] += ( I3 / 2.0 + I4 / du) * ri**2# * r1
            Q[i, j]     += ( I3 / 2.0 - I4 / du) * ri**2# * r2
    
        def f5(rx):
            return Qs(ri, rx)
        def f6(rx):
            return Qs(ri, rx) * rx
    
        # Integrate from 0 to r[0] by linear extrapolation
        r1 = r[0]
        r2 = r[1]
        dr = r2 - r1
        I5, eps5 = integrate.quad(f5, 0.0, r[0])
        I6, eps6 = integrate.quad(f6, 0.0, r[0])
        Q[i, 0] += (  I5*r2/dr - I6/dr) * r1 
        Q[i, 1] += ( -I5*r1/dr + I6/dr) * r2
    
        def f7(rhox):
            return Qs(ri, 1.0/rhox) / rhox**3
        def f8(rhox):
            return Qs(ri, 1.0/rhox) / rhox**2
    
        #
        # Integrate from r[-1] to infinity, 
        # extrapolating U(r) as U[-1] r[-1] / r
        #
        x1 = 1.0/r[-1]
        x2 = 1.0/r[-2]
        I7, eps7 = integrate.quad(f7, 0.0, x1)
        I8, eps8 = integrate.quad(f8, 0.0, x1)
        drho = x2 - x1
        Q[i, -1] += (   + I7*x2/drho - I8/drho)  
        Q[i, -2] += (   - I7*x1/drho + I8/drho)  
        #print i, "**",  Q[i, 0], Q[i, 1], Q[i, -1], Q[i, -2]
        #ff = open("xxx-%d.dat" % i, 'w')
        #for j in range(0, len(r)):
        #    ff.write("%d %g\n" % (j, Q[i, j]))
        #ff.close()
        s0 = np.dot(Q[i, :], 1.0/r)
        s2 = np.dot(Q[i, :], r/r)
        Q[i, i] -= r[i] * s0  
        s = np.dot(Q[i, :], 1.0/r)
        s1 = np.dot(Q[i, :], r/r)
        print "s: ", s, s1, s0, s2
    return Q

def RPA_inter(r):
    """
       Attempt to load interband kernel from file, 
       recalculate if the data is not available
    """
    Rmin = r.min()
    Rmax = r.max()
    N = len(r)
    fname = "rpakernel-inter-Rmin=%g-Rmax=%g-N=%g.dat.npz" % (Rmin, Rmax, N)
    try: 
        data = np.load(fname)
        print "Interband kernel loaded from", fname
        assert linalg.norm(r - data['r']) < 1e-6, "r grids are different"
        return data['Q']
    except:
        import traceback
        traceback.print_exc()
        print "cannot load data from", fname, "; recalculating"
        Q = do_RPA_inter(r)
        np.savez(fname, r=r, Q=Q)
        return Q

if __name__ == '__main__':
   if False:
      testG([1.0, 2.0, 3.0, 5.0, 10.0], np.arange(0.01, 10.0, 0.1), 1.0)
   if False:
      testF(1.0, 2.0, 1.0)
      testF(1.0, 10.0, 1.0)
      testF(1.0, 1.1, 1.0)
   if False:
      show_pi()
   
   rmin = 0.01
   rmax = 100.0
   N = 200
   r = rmin * np.exp(math.log(rmax/rmin)/(N - 1.0) * np.arange(0, N, 1.0))
   #r = np.arange(0.01, 10.0, 0.01)
   kF = 1.0
   Q_inter = RPA_inter(r)
   Q_intra  = RPA_intra(r, kF)
   r_0 = 1.0
   U = 1.0 / (r**2 + r_0**2)**0.5
   def Uq(q):
       return 2.0 * math.pi / (q + 1e-8) * math.exp(-q*r_0)
   
   def rho_direct(x):
       s = 1.0; #max(1.0, kF * x)
       def f(q):
           return Uq(q/s) * Pi_full(q/s, kF) * special.jn(0, q * x / s) * q / s**2
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
   
   rho_TH = -1.0 / 16.0 * r_0 / (r**2 + r_0**2)**1.5
   rho_intra = np.dot(Q_intra, U)
   rho_th_intra = -U * abs(kF) / 2.0 / math.pi
   rho_d = np.vectorize(rho_direct)(r)
   import pylab
   pylab.plot (r, rho_RPA, label='Numerics-kf=0')
   pylab.plot (r, rho_RPA + rho_intra, label='Numerics-kf=%g' % kF)
   pylab.plot (r, rho_TH,  label='Exact-kf=0')
   pylab.plot(r, -U*abs(kF)/2.0/math.pi, label='TF kf=%g' % kF)
   pylab.plot(r, rho_d, label='Direct integration')
   pylab.legend()
   pylab.figure()
   pylab.loglog (r, np.abs(rho_RPA), label='Numerics')
   pylab.loglog (r, np.abs(rho_RPA + rho_intra), label='Numerics-kf=%g' % kF)
   pylab.loglog (r, np.abs(rho_TH),  label='Exact')
   pylab.plot(r, abs(U)*abs(kF)/2.0/math.pi, label='TF kf=%g' % kF)
   pylab.loglog(r, np.abs(rho_d), label='Direct integration')
   pylab.legend()
   pylab.show()
   
