import numpy as np
from math import *
import coulomb

def solve_coulomb(rho_U, U0, r, tau_U, tau_rho):
    M = coulomb_kernel(r)
    U = array(U0)
    rho = np.zeros((len(r)))
    while true:
        U1 = np.dot(M, rho) + U0
        err_U = norm(U - U1)
        U += tau_U * (U1- U)
        rho1 = rho_U(U, r)
        err_rho = norm(rho1 - rho)
        rho += tau_rho * (rho1 - rho)
        print "U error" err_U, "rho error" err_rho
        if (err_U < (1e-06)) and (err_rho < (1e-06)):
            break

    plot(r, U, label="U")


def rho_minus12(U,r):
    return = U / r / 2 / np.pi**2

U0 = (r**2 + r_0**2)**(-0.5)
tau_u = 0.1
tau_rho = 0.1

solve_coulomb(rho_minus12, U0, r, tau_u, tau_rho) 
#does this rho_minus12 need arguements?
