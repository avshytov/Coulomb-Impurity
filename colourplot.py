import numpy as np
import matplotlib
from matplotlib import pylab
from pylab import *

Ustr = 1.3
Bs = []
for Q in range (0,106):
    Bs.extend([0.03**2 * Q**2])
Bs = np.array(Bs)
ivals = [t[0] for t in enumerate(Bs)]
E = np.load('Evec.npy')
r = np.load('rvec.npy')
si = int (1/ r[(len(r)-1)] * len(r))
ldosH = np.zeros((len(E), len(Bs)))
zero = np.zeros((len(Bs)))
first = np.zeros((len(Bs)))
first = np.sqrt(2)*Bs
third = np.zeros((len(Bs)))
third = np.sqrt(3)*Bs
print len(r), 'r', len(E), 'E', len(ldosH), 'ldosH', len(ldosH[0]), 'ldos[0]', len(Bs), 'Bs'

for t in ivals:
    H = Bs[t]
    ldos = np.load('ldosmat-U=%g-m=1-B=%g-grid=%d.npy' %(Ustr, H, len(r)))
    ldosH[:,t] = ldos[:,si]    
print "Producing colour plot..."
np.save('ldos-U=1.3-m=0-N=350',ldosH)
print "data saved"
pcolor(np.sqrt(Bs), E ,ldosH, vmin=0.0, vmax=1.5)
colorbar()
plot(Bs,zero, 'w--')
plot(Bs,first, 'w--')
plot(Bs,-first, 'w--')
plot(Bs,third, 'w--')
plot(Bs,-third, 'w--')
plot(Bs,Bs, 'w--')
plot(Bs,-Bs, 'w--')
ylim(-3.0, 3.0)
xlim(0.0, np.sqrt(Bs[105]))
title('LDOS for U=%g' %Ustr)
xlabel('Square Root of Magnetic Field')
ylabel('Energy')
savefig("ldoscp.pdf")
show()
figure()
