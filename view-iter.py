import glob
import sys
import numpy as np
from pylab import * 
from scipy import linalg

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
   Ndata = min(2, len(datasets))
   last = datasets[-(Ndata):]
   #data = np.load(last[0])
   datas = [] 
   for fname, it in last:
       d = np.load(fname)
       datas.append((it, d))
   #print "available:", data.keys()
   def plot_data(data, key, id, is_normal):
       print data
       if key in data.keys():
          if is_normal:
             plot(data['r'], data[key], label=key + " : " + id)
          else:
             loglog(data['r'], abs(data[key]), label=key + " : " + id)
   def make_figure(datas, keys, is_normal=True):
       figure()
       title ("Z = %g N = %d" % (Z, len(datas[-1][1]['r'])))
       for k in keys:
           for i, d in datas:
               print i
               plot_data(d, k, "t = %d" % i, is_normal)
       legend()
   def draw_slopes(data, keyref):
       r = data['r']
       r_ref = 10.0
       rmin = 3.0
       rmax = max(r)
       i_ref = np.abs(r - r_ref).argmin()
       y_ref = abs(data[keyref][i_ref])
       rvals = [t for t in r if t > rmin and t < rmax]
       xivals = r[i_ref] / rvals
       y1vals = y_ref * xivals
       y2vals = y_ref * xivals**2
       y3vals = y_ref * xivals**3
       loglog(rvals, y1vals, label='1/r')
       loglog(rvals, y3vals, label='1/r^3')
       legend()
   def plot_diff(data, key1, key2):
       Y1 = data[key1]
       Y2 = data[key2]
       dY = Y1 - Y2
       plot(data['r'], dY, label='%s - %s' % (key1, key2))
   make_figure(datas, ['U', 'U1'], False)
   draw_slopes(datas[-1][1], 'U')
   make_figure(datas, ['U', 'U1'], True)
   plot_diff(datas[-1][1], 'U', 'U1')
   legend()
   r = datas[-1][1]['r']
   #U0 = 1.0
   r_0 = 1.0
   Nf = 4
   alpha = 1.0
   eps_th = 1.0 + 3.14159*alpha * Nf / 8
   Uth = Z / np.sqrt(r**2 + r_0**2) / eps_th
   rho_th = -Z * Nf * alpha / eps_th / np.sqrt(r**2 + r_0**2)**3 / 16
   plot(r, Uth - Uth[-1], label='Uth')
   legend()
   make_figure(datas, ['rho', 'rho1'], False)
   draw_slopes(datas[-1][1], 'rho')
   make_figure(datas, ['rho', 'rho1'], True)
   plot(r, rho_th, label='rho_th')
   plot_diff(datas[-1][1], 'rho', 'rho1')
   legend()
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
