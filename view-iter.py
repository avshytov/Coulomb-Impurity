import glob
import sys
import numpy as np
from pylab import * 
from scipy import linalg

if __name__ == '__main__':
   # Search for files matching the template
   fnames = [] 
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
   # Sort data sets
   datasets.sort(lambda x, y: cmp(x[1], y[1]))
   last = datasets[-1]
   data = np.load(last[0])
   print "available:", data.keys()
   def plot_data(key):
       if key in data.keys():
          plot(data[key], label=key)
   def make_figure(keys):
       figure()
       for k in keys:
           plot_data(k)
       legend()
   make_figure(['U', 'U1'])
   make_figure(['rho', 'rho1'])
   iters = [] 
   err_u_data = []
   err_rho_data = []
   for fname, it in datasets:
       print it, fname
       data = np.load(fname)
       U = data['U']
       U1 = data['U1']
       rho = data['rho']
       rho1 = data['rho1']
       err_U = linalg.norm (U - U1)
       err_rho = linalg.norm (rho - rho1)
       iters.append(it)
       err_u_data.append(err_U)
       err_rho_data.append(err_rho)
   figure()
   plot (np.array(iters), np.array(err_u_data), label='err_U')
   plot (np.array(iters), np.array(err_rho_data), label='err_rho')
   legend()
   show()
