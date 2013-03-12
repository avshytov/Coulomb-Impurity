import numpy as np
import matplotlib
from matplotlib import pylab
from pylab import *

Ustr = 1.3
Bs = []
for Q in range (0,106):
    Bs.extend([0.03**2 * Q**2])
ivals = [t[0] for t in enumerate(Bs)]
E = np.load('Evec.npy')
r = np.load('rvec.npy')
si = int (1/ r[(len(r)-1)] * len(r))
ldosH = np.zeros((len(E), len(Bs)))
zero = np.zeros((len(Bs)))
for t in ivals:
    H = Bs[t]
    ldos = np.load('ldosmat-U=%g-m=1-B=%g-grid=200.npy' %(Ustr, H))
    ldosH[:,t] = ldos[:,si]
    
print "Producing colour plot..."
Bs = np.array(Bs)
np.save('ldos-U=1.3-m=0-N=200',ldosH)
print "data saved"
pcolor(np.sqrt(Bs), E ,ldosH, vmin=0.0, vmax=0.5)
colorbar()
plot(Bs,zero, 'w--')
ylim(-3.0, 3.0)
xlim(0.0, np.sqrt(Bs[105]))
title('LDOS for U=%g' %Ustr)
xlabel('Square Root of Magnetic Field')
ylabel('Energy')
savefig("ldoscp.pdf")
show()
figure()
