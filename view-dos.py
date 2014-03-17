import glob
import sys
import numpy as np
from pylab import * 
from scipy import linalg
import odedirac
import util

rc('text', usetex=True)
rc('font', **{'family':'serif', 'serif':['Computer Modern Roman']})

#
#  A script to view iteration progress. Usage:
#
#  python view-iter.py rhoand*npz 
#
#  You can specify a different, more specific template 
#

if __name__ == '__main__':
   # Search for files matching the template
   fnames = [] 
   Z = 0.0
   for arg in sys.argv[1:]:  # there might be more than one arg
       print "arg:", arg
       fnames_arg = glob.glob(arg)
       print "fnames: ", fnames_arg
       fnames.extend(fnames_arg)
   datasets = []
   # Extract iteration number from each file name
   for f in fnames:
       ff = f[:-4] # cut .npz
       fields = ff.split('-')
       print ff, fields
       for l in range(len(fields)):
           field = fields[l]
           print field
           if field[:3] == 'it=':
              iter = int(field[3:])
              datasets.append((f, iter))
           if field[:2] == 'Z=':
              Z = float(field[2:])
           if field[:3] == 'Ef=':
              if len(field) != 3:
                 Ef = float(field[3:])
              else:
                 print 'negative Ef', fields[l+1]
                 Ef = -float(fields[l+1])
   # Sort data sets
   datasets.sort(lambda x, y: cmp(x[1], y[1]))
   
   f = datasets[-1][0]
   data = np.load(f)
   r = data['r']
   U = data['U']

   print fields

   Mmax = 50

   Econv = 1e9 * 1.05457173e-34 * 1e6 / 1.60217657e-19 ###eV                              
   print 'Conversion to eV factor', Econv

   if True:
      r2 = np.linspace(r[0], 10.0, 500)
      U = util.gridswap(r, r2, U)
      r = np.array(r2)
      Espace = np.linspace(-0.8, 0.8, 4000)
      mlist = np.arange(0.0, r[-1]*Espace[-1], 1.0)
      print 'Max m =', mlist[-1]
   else:
      Espace = np.arange(-1.999, 1.999, 0.001) + 1e-5
      mlist = np.arange(0.0, 10.5, 1.0)
   dos = odedirac.doscalc(Espace, r, U, mlist)
   
   Econv = 1e9 * 1.05457173e-34 * 1e6 / 1.60217657e-19 ###eV 

   print 'Conversion to eV factor', Econv
   Espace *= Econv
   dos /= Econv

   figure()
   def _sgn(E):
       if E > 0: return 1.0
       return -1.0; 
   Esgn = np.vectorize (_sgn)(Espace)
   for ri in [0.1, 1.0, 2.0, 3.0, 5.0]:#, 10.0, 15.0]:
       i = np.abs(r - ri).argmin()
       ipeak = abs(dos[i,:]).argmax()
       plot(Espace, dos[i,:], label='r = %g' %r[i])
#       plot(Espace, -dos[i, :]/dos[i,ipeak]*Esgn, label='r = %g' % r[i])
   legend()
   show()
   figure()
   X, Y = meshgrid(Espace[::10], r)
   pcolor(X, Y, dos[:, ::10], vmin=0.0, vmax=1.0)
   colorbar()
   xlim(-0.5, 0.5)
   ylim(r.min(), r.max())
   xlabel('Energy (eV)')
   ylabel('Distance from Impurity (nm)')
   savefig('plots/dosplot-Z=%g-Ef=%g.pdf' %(Z,Ef))
   show()
   
