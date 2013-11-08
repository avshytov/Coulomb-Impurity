import numpy as np
from math import *
from coulomb import *

N = 1000
r = np.zeros((N))
rmin = 0.01
rmax = 1000
for i in range (0,N):
    r[i] = rmin * math.exp(math.log(rmax/rmin)/(N - 1.0) * i)
    #r[i] = 0.1 * (i) + 0.01

r_0 = 1.0

def solve_coulomb(rho_U, U0, r, tau_U, tau_rho):
    M = coulombkernel(r)
    U = array(U0)
    rho = np.zeros((len(r)))
    while True:
        U1 = np.dot(M, rho) + U0
        err_U = norm(U - U1)
        U += tau_U * (U1- U)
        rho1 = rho_U(U, r)
        err_rho = norm(rho1 - rho)
        rho += tau_rho * (rho1 - rho)
        print "U error", err_U, "rho error", err_rho
        if (err_U < (1e-10)) and (err_rho < (1e-10)):
            break

    return U

Nf = 4
alpha = 2.5

def rho_minus12(U,r):
    return -U / r / 2.0 / np.pi**2 * Nf * alpha

U0 = (r**2 + r_0**2)**(-0.5)
tau_u = 0.1
tau_rho = 0.1

U = solve_coulomb(rho_minus12, U0, r, tau_u, tau_rho) 

ra = 10.0
rb = 33.0

ivals = [t[0] for t in enumerate(r) if t[1] < rb and t[1] > ra]
xvals = [r[t] for t in ivals]
yvals = [(abs(U[t])) for t in ivals]

grad = np.polyfit(log(xvals), log(yvals), 1)
print grad
fitvals = exp(grad[1] + grad[0]*log(xvals))

figure()
plot (r, U)
figure()
loglog (r, abs(U), label='U')
loglog (r, 1.0/r, label='1/r')
loglog(r, 1.0/r/r, label='1/r^2')
loglog(xvals, abs(fitvals), label='Fit: alpha = %g' % (-grad[0]))
legend()
show()
