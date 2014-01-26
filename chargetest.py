import density
from density import *
import pylab as pl

rmin = 0.01
rmax = 50.0
N = 500
r = util.make_lin_grid(rmin, rmax, N)
Ef = 0.3
Ecut = -3.0
graphene = GrapheneResponse(r, Ef, Ecut=Ecut, Mmax=31)

Emax = 2.0
Es = np.arange(-Emax, Emax, 0.01)

Z = 0.25
r_0 = 1.0
U = Z / np.sqrt(r**2 + r_0**2)
rho_th = -Z * r_0 / 16.0 / np.sqrt(r**2 + r_0**2)**3

#rho_RPA = graphene.rho_RPA(U)
#rho_U   = graphene.rho_U(U)
#rho_Ub  = graphene.bandResponse(U)
#rho_RPAb = graphene.apply_kernel(graphene.Q_Emax - graphene.Q_Emin, U)
rho_0 = (graphene.Emax**2 - graphene.Emin**2) / 4.0 / math.pi
#rho_1 = graphene.diracDensity(U)
mlist = graphene.mlist
imin = 0
imax = len(r) - 1
rmin = r[0]
rmax = r[-1]
dr = r[imin+1]-r[imin]
Ugrid = []
Qtots = []
Qths = []

print 'Total Charge Calculation: rmin=', rmin, 'rmax=', rmax
for Z0 in [0.02, 0.1, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4,
           1.8, 2.2, 2.6, 3.0]: 

    print 'Calculating for Z0=',Z0
    Ugrid.append(Z0)
    theorytest = -Z0 / 16.0 / (r**2 + r_0**2)**(1.5)
    Qtheory = (( 1.0 / np.sqrt(r_0**2+rmin**2))
                - (1.0/ np.sqrt(r_0**2 + rmax**2)))
    Qtheory *= Z0 * np.pi / 8.0 * r_0
    Us = Z0 / np.sqrt(r**2 + r_0**2)
    #rhotot = graphene.rho_U(Us)
 
    
    rho_RPAbt = graphene.apply_kernel(graphene.Q_Emax - graphene.Q_Emin, Us)
    rho_Ub  = graphene.bandResponse(Us)
    rho_rpa_max = graphene.rho_RPA(Us)
    rho_rpa_min = graphene.apply_kernel(graphene.Q_Emin, Us)
    rhotot = rho_Ub + rho_rpa_min
    
    theorytest = rho_rpa_max
    
#    pl.figure() 
#    pl.title('Total density for Z0=%g'%Z0)
#    pl.plot(r,rhotot, label='total charge density')
#    pl.plot(r,theorytest, '--', label='theory test')
#    pl.plot(r,theorytest - rhotot, '--', label='diff')
#    pl.legend()

#    pl.figure()
#    pl.title('loglog Total density for Z0=%g'%Z0)
#    pl.loglog(r,abs(rhotot/rhotot[0]), label='total charge density')
#    pl.loglog(r,abs(theorytest), label='theory test')
#    pl.legend()

#    pl.figure() 
#    pl.title('Total band density for Z0=%g'%Z0)
#    pl.plot(r,rho_Ub, label='total charge density')
#    pl.plot(r,rho_RPAbt, '--', label='theory test')
#    pl.plot(r,rho_RPAbt - rho_Ub, '--', label='diff')
#    pl.legend()

    #pl.figure()
    pl.title('loglog Total band density')
    pl.loglog(r,abs(rho_Ub/rho_Ub[0]), label='total charge density Z0=%g'%Z0)
    pl.loglog(r,abs(rho_RPAbt/rho_RPAbt[0]), '--')
#    pl.loglog(r,1/r**2, '--', label='1/r^2')
#    pl.loglog(r, 1/r**3, '--', label='1/r^3')
    pl.legend()

    Qsim = 0.5 * dr * rhotot[imin] * r[imin]
    Qsim += 0.5 * dr * rhotot[imax] * r[imax]
    for i in range (imin+1, imax):
        Qsim += dr * rhotot[i] * r[i]
    Qsim *= -2.0 * np.pi
    Qtots.append(Qsim)
    Qths.append(Qtheory)

Ugrid = np.array(Ugrid)
Q_sim = np.array(Qtots)
Q_linear = np.array(Qths)
Q_bs = (np.pi/8.0*Ugrid + 0.19*Ugrid**3)
pl.loglog(r,1/r**2, '--', label='1/r^2')
pl.loglog(r,1/r**3, '--', label='1/r^3')

#pl.figure()
#pl.plot(Ugrid, Q_sim, label='sim')
#pl.plot(Ugrid, Q_linear, label='Linear')
#pl.plot(Ugrid, Q_bs, label='B-S')
#pl.legend()
pl.show()
