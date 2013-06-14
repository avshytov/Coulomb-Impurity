import numpy as np
import matplotlib
from matplotlib import pylab
from pylab import *

Ustr = 1.3
Bs = np.load('Bs-106.npy')
#Bs=[2.3409,2.34098,2.34122,2.34162,2.34218,2.3429,2.34378,2.344,2.3444,2.34482,2.34602,2.34738,2.3489,2.35058,2.35242,2.35442,2.35658,2.3589,2.36138,2.36402,2.36682,2.3696,2.36978,2.37]
Bs = np.array(Bs)
ivals = [t[0] for t in enumerate(Bs)]
E = np.load('Evec-N=200.npy')
r = np.load('rvec-N=200.npy')
info = np.load('EMinfo.npy')
si = int (3.0/ r[(len(r)-1)] * len(r))
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
np.save('ldos-U=%g-m=0-N=%d-Bs=%d' %(info[0],len(r),len(Bs)),ldosH)
print "data saved"
pcolor(np.sqrt(Bs), E ,ldosH, cmap=cm.hot, vmin=0.0, vmax=0.50)
colorbar()
plot(Bs,zero, 'b--')
plot(Bs,first, 'g--')
plot(Bs,-first, 'g--')
#plot(Bs,third, 'w--')
#plot(Bs,-third, 'w--')
#plot(Bs,Bs, 'w--')
#plot(Bs,-Bs, 'w--')
ylim(-6.0, 6.0)
xlim(np.sqrt(Bs[0]), np.sqrt(Bs[len(Bs)-1]))
title('LDOS-H for U=%g, r=%g' %(Ustr, r[si]))
xlabel('Square Root of Magnetic Field')
ylabel('Energy')
savefig("ldoscp.pdf")
show()
figure()
