import numpy as np
from numpy import linalg
import math
import scipy
import mkrpa2
import util
import iteration
import density
import pylab as pl
import coulomb

N = 500
Mmax = 31
rmin = 0.01
rmax = 50.0

r = util.make_lin_grid(rmin, rmax, N)
Ef = 0.3
B0 = 0.0
Ecut = -3.0
Z = 3
r_0 = 1.0
Nf = 4
alpha = 0.9


tau_u_set = iteration.make_tau_set(3, 0.2, 0.1)
print "taus: ", tau_u_set
tau_rho_set = list(tau_u_set)
Uext = -Z * alpha / np.sqrt(r**2 + r_0**2)
rho_ext = -Z * alpha / 16.0 / np.sqrt(r**2 + r_0**2)**3  

graphene = density.GrapheneResponse (r, Ef, Ecut=Ecut, Mmax=Mmax)


def rho_U(U, r):
    return graphene.rho_U(U) * Nf * alpha

def rho_RPA(U, r):
    rho = graphene.rho_RPA(U) * Nf * alpha
    #rho[-50:-1] = 0.0 
    return rho

#rho_rpa = graphene.rho_RPA(Uext)
#C = coulomb.kernel(graphene.rexp)
#Unew = graphene.apply_kernel(C, rho_rpa) * 8.0 / math.pi * (-1)
#pl.plot (r, Uext, label='Uext')
#pl.plot (r, rho_RPA(Uext, r), label='rpa')
#pl.plot (r, rho_ext, label='rpa exp')
#pl.plot (r, Unew, label='Unew')
#pl.legend()
#pl.show()

fname_templ = util.makeFileTemplate("data/denssim-Z=%g-N=%g-Nfa=%g-Ef=%g" 
                                    % (Z, N, Nf*alpha, Ef))

iteration.solve_coulomb(rho_U, Uext, r, tau_u_set, tau_rho_set, 
                        fname_template=fname_templ)

