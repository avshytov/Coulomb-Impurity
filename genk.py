#
#
import rpakernel
import util
#import pylab as pl
import numpy as np
import mkrpa2
import rpam
import rpacorr

#kF = 3.0
rmin = 0.01
rmax = 50.0
N = 500
Mmax = 31
mlist = np.array(range(0, Mmax))
rexp = util.make_exp_grid(rmin, rmax, N)
rlin = util.make_lin_grid(rmin, rmax, N)

#Q0 = rpakernel.kernel_m_inter(rexp, mlist, 'exp')
#Q1 = mkrpa2.RPA_inter(rexp, 'exp')
Q2 = rpacorr.RPA_corr_inter(rexp, 'exp')

#for kF in [3.0, 2.0, 1.0, 0.5, 0.3, 0.2, 0.1, 2.5, 1.5, 0.5, 0.4, 0.35, 0.25, 0.15, 0.05]:
#for kF in [3.0, 1.0, 0.5, 0.4, 0.3, 0.2, 0.1]:
#for kF in np.arange(0.01, 0.5, 0.01):
grid_intra = 'exp'

r_dict = {
  'lin' : rlin, 
  'exp' : rexp
}

intra_r = r_dict[grid_intra]

for kF in np.arange(0.5, 0.01, -0.01):
#for kF in [1.0, 0.3, 0.1, 2.5, 1.5]:
#for kF in [2.5, 1.5]:
#for kF in [0.1, 0.2, 0.3, 0.4, 0.05, 0.15, 0.25, 0.35, 0.5]:
#for kF in [3.0, 0.1, 0.2, 0.3]:
#for kF in [1.0, 0.05, 0.15, 0.25, 0.35]:
#for kF in [1.0]:
#for kF in [0.05, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.4, 0.5]:
    #Q_tot = mkrpa2.RPA_intra(intra_r, kF, grid_intra)
    #Q_m   = rpakernel.kernel_m_intra(rlin, mlist, kF, '4')
    #Q_m   = rpam.kernel_m_intra(intra_r, mlist, kF, '3' + grid_intra)
    Q_c   = rpacorr.RPA_corr_intra(intra_r, kF, grid_intra)
