#
#
import rpakernel
import util
#import pylab as pl
import numpy as np
import mkrpa2
import rpam

#kF = 3.0
rmin = 0.01
rmax = 50.0
N = 500
Mmax = 21
mlist = np.array(range(0, Mmax))
rexp = util.make_exp_grid(rmin, rmax, N)
rlin = util.make_lin_grid(rmin, rmax, N)

Q0 = rpakernel.kernel_m_inter(rexp, mlist)
Q1 = mkrpa2.RPA_inter(rexp)

#for kF in [3.0, 2.0, 1.0, 0.5, 0.3, 0.2, 0.1, 2.5, 1.5, 0.5, 0.4, 0.35, 0.25, 0.15, 0.05]:
#for kF in [1.0, 0.3, 0.1, 2.5, 1.5]:
#for kF in [2.5, 1.5]:
#for kF in [0.1, 0.2, 0.3, 0.4, 0.05, 0.15, 0.25, 0.35, 0.5]:
#for kF in [3.0, 0.1, 0.2, 0.3]:
#for kF in [1.0, 0.05, 0.15, 0.25, 0.35]:
for kF in [3.0]:
    Q_tot = mkrpa2.RPA_intra(rlin, kF)
    #Q_m   = rpakernel.kernel_m_intra(rlin, mlist, kF)
    Q_m   = rpam.kernel_m_intra(rlin, mlist, kF)
    
