import rpam
import rpakernel

import numpy as np
r = np.linspace(0.01, 50.0, 500)
m = np.arange(0.0, 31.0, 1.0)
print r
print m
kF = 3.0

Q = rpam.kernel_m_intra(r, m, kF, '3')

#Q1 = rpakernel.kernel_m_intra(r, m, kF)
Q1 = rpam.kernel_m_intra(r, m, kF, '4')

import pylab as pl

for r_i in [0.1, 1.0, 2.0, 5.0, 7.0]:
    i = np.abs(r - r_i).argmin()
    pl.figure()
    pl.plot(r, Q[i, :], label='Q0: r = %g' % r[i])
    pl.plot(r, Q1[i, :], label='Q1: r = %g' % r[i])
    pl.plot(r, Q1[i, :] - Q[i, :], label='diff: r = %g' % r[i])
    pl.legend()
pl.show()

