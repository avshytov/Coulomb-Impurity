import numpy as np
import matplotlib
from matplotlib import pylab
from pylab import *

E = np.load('Evec.npy')
r = np.load('rvec.npy')
a0 = np.load('ldosmat-U=0-m=0.npy')
a1 = np.load('ldosmat-U=0-m=1.npy')
a2 = np.load('ldosmat-U=0-m=2.npy')
a3 = np.load('ldosmat-U=0-m=3.npy')
b0 = np.load('ldosmat-U=sub-m=0.npy')
b1 = np.load('ldosmat-U=sub-m=1.npy')
b2 = np.load('ldosmat-U=sub-m=2.npy')
b3 = np.load('ldosmat-U=sub-m=3.npy')
c0 = np.load('ldosmat-U=sup-m=0.npy')
c1 = np.load('ldosmat-U=sup-m=1.npy')
c2 = np.load('ldosmat-U=sup-m=2.npy')
c3 = np.load('ldosmat-U=sup-m=3.npy')


for rs in [1.0]:
            si = int (rs / r[(len(r)-1)] * len(r))
            ri = r[si]
            
            print "plotting zero potential LDOS"
            plot (E, a0[:, si], label='LDOS, m=0, r = %g' % ri)
            plot (E, a1[:, si], label='LDOS, m=1, r = %g' % ri)
            plot (E, a2[:, si], label='LDOS, m=2, r = %g' % ri)
            plot (E, a3[:, si], label='LDOS, m=3, r = %g' % ri)
            title('With no external potential')
            legend()
            xlim(-2.5, 2.5)
            figure()

            print "plotting supercritical LDOS"
            plot (E, b0[:, si], label='LDOS, m=0, r = %g' % ri)
            plot (E, b1[:, si], label='LDOS, m=1, r = %g' % ri)
            plot (E, b2[:, si], label='LDOS, m=2, r = %g' % ri)
            plot (E, b3[:, si], label='LDOS, m=3, r = %g' % ri)
            title('With subcritical potential')
            legend()
            xlim(-2.5, 2.5)
            figure()
            

            print "plotting supercritical LDOS" 
            plot (E, c0[:, si], label='LDOS, m=0, r = %g' % ri)
            plot (E, c1[:, si], label='LDOS, m=1, r = %g' % ri)
            plot (E, c2[:, si], label='LDOS, m=2, r = %g' % ri)
            plot (E, c3[:, si], label='LDOS, m=3, r = %g' % ri)
            title('With supercritical potential')
            legend()
            xlim(-2.5, 2.5)
show()
