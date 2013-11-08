import numpy as np
import math

def propfit(y, x):
    sxy = np.sum(y * x)
    sx2 = np.sum(x**2)
    a = sxy / sx2
    return a


def linfit(y, x):
    n = len(x)
    sx = np.sum(x) / n
    sy = np.sum(y) / n
    sx2 = np.sum(x**2) / n
    sxy = np.sum(x * y) / n

    a = (sxy - sx * sy) / (sx2 - sx**2)
    b = (sy - a * sx);

    return a, b


if __name__ == '__main__': 
   x = np.array([1, 2, 3, 4, 5])
   import random
   err = x * 0.0
   for i in range(len(x)):
       err[i] = 0.9 * (2 * random.random() - 1.0)
   y = 2.1 * x + err
   
   print "x = ", x
   print "y = ", y
   a_prop = propfit(y, x)
   print "propfit: ", a_prop
   from pylab import figure, plot, legend, show, title
   plot(x, y, 'ko', label='data')
   plot(x, a_prop * x, label='fit')
   legend()
   title('Proportional fit')
   #show()

   figure()
   y = 2.7 * x - 1.1 + err
   print "x = ", x
   print "y = ", y
   a, b = linfit(y, x)
   print "linfit: ", a, b
   from pylab import figure, plot, legend, show
   plot(x, y, 'ko', label='data')
   plot(x, a * x + b, label='fit')
   legend()
   title('Linear fit')
   show()
    
