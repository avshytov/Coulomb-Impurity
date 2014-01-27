#from numpy import * 
import math
import RPA
from scipy import integrate
import numpy as np
from scipy import linalg
from scipy import special

integrate_all = False

def do_kernel_m_intra_new(r, mvals, kF):
    N = len(r)
    r_out = r
    r_in = []
    dr0 = 0.1 / kF
    Q1 = np.zeros((N, N))
    for i in range(N - 1):
        r_in.append(r[i])
        dr = r[i + 1] - r[i]
        nsplit = int( dr / dr0 ) + 1
        for k in range(1, nsplit):
            r_in.append(r[i] + dr/nsplit * k)
    #Nk = 100 * kF * r.max()
    r_in = np.array(r_in)
    def F(r1, r2):
        ra, rb = r1, r2
        if (r2 < r1):
            ra, rb = r2, r1
        def G(k):
            s = 0.0
            for m in mvals:
                j1a = special.jn(m,     k * ra)
                j1b = special.jn(m,     k * rb)
                j2a = special.jn(m + 1, k * ra)
                j2b = special.jn(m + 1, k * rb)
                y1b = special.yn(m,     k * rb)
                y2b = special.yn(m + 1, k * rb)
                s += (j1a**2 + j2a**2) * (j1b * y1b + j2b * y2b)
            s *= k**2
            return s
        if True:
           N0 = 50
           Nk = 2 * int(((r1 + r2) * kF + 1) * N0) + 1
           #print r1, r2, Nk
           ks = np.linspace(1e-4/r.max(), kF, Nk)
           Gs = G(ks)
           weights = np.zeros(np.shape(ks))
           weights[0::2] = 2.0
           weights[1::2] = 4.0
           weights[0] = weights[-1] = 1.0
           weights /= 3.0
           dk = ks[1] - ks[0]
           I = np.dot(weights ,  Gs) * dk
        #return np.dot(Gs, weights) / 8.0 / math.pi * dk
        else:
           I, eps = integrate.quad(G, 0.0, kF)
           print r1, r2
        return I / 8.0 / math.pi
    for i in range(N):
        ri = r[i]
        #def Gfun(r):
        #    return Qsum(ri, r2)
        for j in range(N - 1):
            r1 = r[j]
            r2 = r[j + 1]
            rc = (r1 + r2) / 2.0
            dr = r2 - r1
            nj = int(dr / dr0) + 1
            Gk = np.zeros((nj,))
            xk = np.linspace(r1, r2, nj + 2)[1:-1]
            for k in range(nj):
                Gk[k] = F(r[i], xk[k]) * xk[k]
            dxk = dr / nj
            I1 = sum(Gk) * dxk
            I2 = sum(Gk * (xk - rc)) * dxk
            Q1[i, j    ] += 0.5 * I1 - I2 / dr
            Q1[i, j + 1] += 0.5 * I1 + I2 / dr
        Q1[i, 0] += F(r[i], 0) * r[0]**2 / 2.0
        print "m-intra-new; i = ", i, sum(Q1[i, :]) * 8.0 * math.pi**2
    return Q1 * 4.0 * math.pi 
            
def do_kernel_m_intra(r, mvals, kF):
    N = len(r)
    Q1 = np.zeros((N, N))
    dr = np.zeros((N))
    dr[1:] = np.diff(r)
    dr[0] = r[1] - r[0]
    dr0 = min(0.5 / kF, 0.5)
    def Qsum(r1, r2):
        s = 0.0
        #print r1, r2
        for m in mvals:
            s += RPA.Qm1(kF, m, r1, r2)
        return s
    for i in range(N):
        ri = r[i]
        #def Gfun(r):
        #    return Qsum(ri, r2)
        for j in range(N - 1):
            r1 = r[j]
            r2 = r[j + 1]
            rc = (r1 + r2) / 2.0
            dr = r2 - r1
            nj = int(dr / dr0) + 1
            Gk = np.zeros((nj,))
            xk = np.linspace(r1, r2, nj + 2)[1:-1]
            for k in range(nj):
                Gk[k] = Qsum(r[i], xk[k]) * xk[k]
            dxk = dr / nj
            I1 = sum(Gk) * dxk
            I2 = sum(Gk * (xk - rc)) * dxk
            Q1[i, j    ] += 0.5 * I1 - I2 / dr
            Q1[i, j + 1] += 0.5 * I1 + I2 / dr
        Q1[i, 0] += Qsum(r[i], 0) * r[0]**2 / 2.0
        print "m-intra; i = ", i, sum(Q1[i, :]) * 8.0 * math.pi**2
    return Q1 * 4.0 * math.pi 
            
def do_kernel_m_inter(r, mvals):
    Qs = RPA.Qs_spline(mvals)
    N = len(r)
    Q = np.zeros ((N, N))
    for i in range (0, N):
        print "Q: ", i
        ri = r[i]
        for j in range (0, i):
            r1 = r[j]
            r2 = r[j + 1]
            dr = r2 - r1
            rc = (r1 + r2)/2.0
            def f1(rx):
                return Qs(ri, rx)
            def f2(rx):
                return Qs(ri, rx) * (rx - rc)
            if integrate_all or (abs(i - j) < 10) or (i < 20):
                 I1, eps1 = integrate.quad(f1, r1, r2)
                 I2, eps2 = integrate.quad(f2, r1, r2)
            else:
                 I1 = Qs(ri, rc) * dr
                 qs1 = Qs(ri, r1)
                 qs2 = Qs(ri, r2)
                 I2 = (qs2 - qs1) * dr**2 / 12.0
            Q[i, j]     += ( I1 / 2.0 - I2 / dr ) * r1 
            Q[i, j + 1] += ( I1 / 2.0 + I2 / dr ) * r2
        
        for j in range (i + 1, N):
            x1 = 1.0/r[j]
            x2 = 1.0/r[j - 1]
            drho = x2 - x1
            r1 = r[j - 1]
            r2 = r[j]
            rhoc = (x1 + x2)/2.0
            def f3(rhox):
                return Qs(ri, 1.0/rhox) / rhox**3 
            def f4(rhox):
                return Qs(ri, 1.0/rhox) / rhox**3 * (rhox - rhoc)
            if integrate_all or (abs(i - j) < 10) or (i < 20):
                I3, eps3 = integrate.quad(f3, x1, x2)
                I4, eps4 = integrate.quad(f4, x1, x2)
            else:
                I3  = Qs(ri, 1.0/rhoc) / rhoc**3 * drho
                qs1 = Qs(ri, 1.0/x1) / x1**3
                qs2 = Qs(ri, 1.0/x2) / x2**3
                I4 = (qs2 - qs1)*drho**2/12.0
            Q[i, j - 1] += ( I3 / 2.0 + I4 / drho) 
            Q[i, j]     += ( I3 / 2.0 - I4 / drho) 
    
        def f5(rx):
            return Qs(ri, rx)
        def f6(rx):
            return Qs(ri, rx) * rx
    
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
    
        x1 = 1.0/r[-1]
        x2 = 1.0/r[-2]
        I7, eps7 = integrate.quad(f7, 0.0, x1)
        I8, eps8 = integrate.quad(f8, 0.0, x1)
        drho = x2 - x1
        #Q[i, -1] += (   + I7*x2/drho - I8/drho)  
        #Q[i, -2] += (   - I7*x1/drho + I8/drho)  
        C1 = (   + I7*x2/drho - I8/drho)  
        C2 = (   - I7*x1/drho + I8/drho)  
        Q[i, -1] += C1#(   + I7*x2/drho - I8/drho)  
        Q[i, -2] += C2# (   - I7*x1/drho + I8/drho)  
        
        s = np.dot(Q[i, :], 1.0/r)
        Q[i, i] -= r[i] * s
        s = np.dot(Q[i, :], 1.0/r)
        s1 = np.dot(Q[i, :], r/r)
        print "s: ", s, s1
        Q[i, -1] -= C1
        Q[i, -2] -= C2
    return Q

def kernel_m_inter(r, mlist):
    fname = "data/rpa-m-Rmin=%g-Rmax=%g-N=%g-Mmax=%d.npz" % (r.min(), r.max(), len(r), mlist[-1])
    try:
        data = np.load(fname)
        print "Loading RPA m-resolved kernel from", fname
        print "Check r"
        assert linalg.norm(r - data['r'])<1e-8, "r vectors match"
        assert linalg.norm(np.array(mlist) - data['mlist'])<1e-8, "m lists match"
        return data['Qm'] / math.pi**2
    except:
        import traceback
        traceback.print_exc()
        print "cannot load", fname, ": recalculating"
        Qm = do_kernel_m_inter(r, mlist) 
        np.savez(fname, Qm=Qm, r=r, mlist=np.array(mlist))
        return Qm / math.pi**2
  
def kernel_m_intra(r, mlist, kF):
    fname = "data/rpa-m2-Rmin=%g-Rmax=%g-N=%g-Mmax=%d-kF=%g.npz" % (r.min(), r.max(), len(r), mlist[-1], kF)
    try:
        data = np.load(fname)
        print "Loading RPA m-resolved intraband kernel from", fname
        print "Check r"
        assert linalg.norm(r - data['r'])<1e-8, "r vectors match"
        assert linalg.norm(np.array(mlist) - data['mlist'])<1e-8, "m lists match"
        return data['Q'] 
    except:
        import traceback
        traceback.print_exc()
        print "cannot load", fname, ": recalculating"
        Q = do_kernel_m_intra(r, mlist, kF) 
        np.savez(fname, Q=Q, r=r, mlist=np.array(mlist))
        return Q
  
def kernel_m (r, mlist, kF):
    Q1 = kernel_m_inter(r, mlist)
    if (kF * r.max() > 0.01):
        Q1 += kernel_m_intra(r, mlist, kF)
    return Q1
