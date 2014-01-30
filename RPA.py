from scipy import special
from scipy import integrate
from pylab import * 
from scipy import interpolate


def Fm(x, m, r1, r2, semi = True):
    ra = r1
    rb = r2
    if (r2 < r1):
       rb = r1
       ra = r2
    if m > 20 and semi:
       
       za = x * ra
       zb = x * rb
       if (za < 1e-5) and (zb < 1e-5): 
           e = -1 + m * math.log(ra/rb * (m + 1.0)/m) -0.5 *math.log((m+1.0)*m) - math.log(rb/2/(m + 1.0))
           return -math.exp(2*e)/4.0
           #return 0.0
       sqa1 = math.sqrt(m**2 + za**2)
       sqa2 = math.sqrt((m + 1)**2 + za**2)
       sqb1 = math.sqrt(m**2 + zb**2)
       sqb2 = math.sqrt((m + 1)**2 + zb**2)
       
       li1 = sqa1 + m * math.log(za/(m + sqa1)) - 0.5 * math.log(sqa1) 
       li2 = sqa2 + (m + 1) * math.log(za/(m + 1 + sqa2)) - 0.5 * math.log(sqa2) 
       lk1 = -sqb1 - m * math.log(zb/(m + sqb1)) - 0.5 * math.log(sqb1)
       lk2 = -sqb2 - (m + 1)* math.log(zb/(m + 1 + sqb2)) - 0.5 * math.log(sqb2)
       L =  2.0*(li1 + lk1)
       Li = 2.0*(li2 - li1)
       Lk = 2.0*(lk2 - lk1)
       return math.exp(L) * x * x * (math.exp(Li) - 1) * (math.exp(Lk) - 1)/4
       #return math.exp(L) * (math.exp(Li) - 1) * (math.exp(Lk) - 1)/4
       
    im1 = special.ive(m, x * ra)
    km1 = special.kve(m, x * rb)
    
    im2 = special.ive(m + 1, x * ra)
    km2 = special.kve(m + 1, x * rb)
    
    e = math.exp(x * (ra - rb))
    
    q1 = im1 * km1 * im1 * km1 
    q2 = im2 * km2 * im2 * km2 
    q3 = im1 * km2 * im1 * km2
    #q4 = q3
    q4 = im2 * km1 * im2 * km1
    #q = x * (im1 * km1 - im2 * km2) * e
    #print r1, r2, m, e, q
    return x * x * (im1 * im1 - im2 * im2) * (km2 * km2 - km1 * km1) * e * e

def Gm(E, m, r1, r2):
    ra = r1
    rb = r2
    if (r2 < r1):
        ra, rb = rb, ra
    j1a = special.jn(m,     E * ra)
    j2a = special.jn(m + 1, E * ra)
    j1b = special.jn(m,     E * rb)
    j2b = special.jn(m + 1, E * rb)
    y1b = special.yn(m,     E * rb)
    y2b = special.yn(m + 1, E * rb)
    return E * E * (j1a * j1a + j2a * j2a) * (j1b * y1b + j2b * y2b) 

if False:
  #for m in [0, 1, 2, 5 10]: # 20, 40]:
  for m in [19, 20, 21, 22]:
    figure()
    for r2x in [1.001, 1.01, 1.1, 2.0]:
        xvals = arange (0.0, 100.0, 0.01)
        S = 1.0
        yvals = vectorize(lambda x: Qm(x, m, 1.0, r2x)*S)(xvals)
        plot (xvals, yvals, label='m = %d r2 = %g' % (m, r2x))
    legend()    
  show ()
  
if False:
  #for m in [0, 1, 2, 5, 10, 20, 40]:
  for m in [21, 22, 30, 40, 50]:
    figure()
    for r2x in [1.001, 1.01, 1.1, 2.0]:
        xvals = arange (0.0, 100.0, 0.01)
        S = 1.0
        yvals1 = vectorize(lambda x: Qm(x, m, 1.0, r2x, False)*S)(xvals)
        yvals2 = vectorize(lambda x: Qm(x, m, 1.0, r2x, True)*S)(xvals)
        plot (xvals, yvals2/yvals1-1, label='m = %d r = %g' % (m, r2x))
        #plot (xvals, yvals1, '--', label='m = %d r2 = %g, true' % (m, r2x))
        #plot (xvals, yvals2, label='m = %d r2 = %g, approx' % (m, r2x))
    legend()    
  show ()

def Qm(m, r1, r2):
    def F(x):
        return Fm(x, m, r1, r2) 
    I, eps = integrate.quad(F, 0.0, Inf, limit=1000)
    #print '**', m, r1, r2, I, eps
   
    return I 

def Qm1(E_F, m, r1, r2):
    def G(E):
        return Gm(E, m, r1, r2)
    Emax = abs(E_F)
    I, eps = integrate.quad(G, 0.0, Emax)
    return I / 8.0 / math.pi 

if False:
    #mvals = range (0, 50)
    rvals = arange(0.01, 3.0, 0.01)
    figure()
    for m in [0, 1, 2, 5]:
    #for r2x in [1.001, 1.01, 1.1, 1.5]:
        #for m in mvals: 
        pmvals = vectorize(lambda r: Qm(m, r, 1.0)*2.0*math.pi)(rvals)
        plot (rvals, pmvals, label=r'$m = %d$' % m)
    xlabel(r"Ratio of the two distances, $r/r'$")
    ylabel(r"Screening kernel, $Q_m(r, r')$")
    legend()
    savefig("qm-1.jpg")
    show()
    

def Q_reg(r1, r2, eps):
    def f(theta):
        R = math.sqrt(eps**2 + r1*r1 + r2*r2 - 2 * r1 * r2 * math.cos(theta))
        return (1.0 - 3.0*eps**2/R**2) / R**3
    theta0 = 5 * eps / (r1 + r2) / 2.0
    if (theta0 > 0.1):
        theta0 = 0.1
    I, epsI = integrate.quad(f, 0.0, math.pi)
    #I1, eps1 = integrate.quad(f, 0.0, theta0)
    #I2, eps2 = integrate.quad(f, theta0, math.pi)
    #I = I1 + I2
    return I / math.pi / 32.0 / math.pi



def Qs_spline(mvals):
    xvals = arange (0.0000, 1.0001, 0.001)
    def Qmsum(x):
        s = 0.0
        for m in mvals:
            s += Qm(m, x, 1.0)
        return s    
    yvals = vectorize(Qmsum)(xvals)
    #import pylab as pl 
    #pl.plot(xvals, yvals)
    #pl.title('m = %g %g' % (mvals[0], mvals[-1]))
    #pl.show()
    spl = interpolate.splrep(xvals, yvals)
    def Qm_s (r1, r2):
        if (r1 <= r2):
            x = r1 / r2
            r = r2
        else:
            x = r2 / r1
            r = r1
        #if (x < x0): 
        #    return 0.0
        y = interpolate.splev (x, spl, der=0)
        return y / r**3  
    
    return Qm_s
   

def Qm_spline(m):
    return Qs_spline([m])
    xvals = arange (0.0001, 1.0001, 0.001)
    yvals = vectorize(lambda x: Qm(m, x, 1.0))(xvals)
    spl = interpolate.splrep(xvals, yvals)
    def Qm_s (r1, r2):
        if (r1 <= r2):
            x = r1 / r2
            r = r2
        else:
            x = r2 / r1
            r = r1
        #if (x < x0): 
        #    return 0.0
        y = interpolate.splev (x, spl, der=0)
        return y / r**3 
    
    return Qm_s
