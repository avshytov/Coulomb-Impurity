import glob
import sys
import numpy as np
from pylab import * 
from scipy import linalg
import odedirac

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
       for field in fields:
           print field
           if field[:3] == 'it=':
              iter = int(field[3:])
              datasets.append((f, iter))
           if field[:2] == 'Z=':
              Z = float(field[2:])
   # Sort data sets
   datasets.sort(lambda x, y: cmp(x[1], y[1]))
   
   f = datasets[-1][0]
   data = np.load(f)
   r = data['r']
   U = data['U']

   Mmax = 50
   Espace = np.arange(-1.999, 1.999, 0.001) + 1e-5
   mlist = np.arange(0.0, 10.5, 1.0)
   dos = odedirac.doscalc(Espace, r, U, mlist)
   figure()
   for ri in [0.1, 1.0, 2.0, 3.0, 5.0, 10.0]:
       i = np.abs(r - ri).argmin()
       plot(Espace, dos[i, :], label='r = %g' % r[i])
   legend()
   figure()
   X, Y = meshgrid(Espace[::10], r)
   pcolor(X, Y, dos[:, ::10])
   colorbar()
   
   show()
   
