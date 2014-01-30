using_c_code = False

if using_c_code: 
   import _rpam
   def do_kernel_m_intra(r, mlist, kF):
       return _rpam.kernel_m_intra(list(r), list(mlist), kF)
                    
else:
   import rpakernel
   def do_kernel_m_intra(r, mlist, kF):
       return rpakernel.do_kernel_m_intra(r, mlist, kF)


import _rpam
import numpy as np
from scipy import linalg

def kernel_m_intra(r, mlist, kF, version='3'):
    fname =  "data/rpa-m%s-Rmin=%g-Rmax=%g-N=%g-Mmax=%d-kF=%g.npz" % (version, r.min(), r.max(), len(r), mlist[-1], kF)
    try: 
        data = np.load(fname); 
        print "Loading RPA m-resolved intraband kernel from", fname
        print "Check r"
        assert linalg.norm(r - data['r'])<1e-8, "r vectors match"
        assert linalg.norm(np.array(mlist) - data['mlist'])<1e-8, "m lists match"
        return data['Q']
    except: 
        import traceback
        traceback.print_exc()
        print "cannot load", fname, ": recalculating"
        Q = _rpam.kernel_m_intra(list(r), list(mlist), kF)
        Q = np.array(Q)
        np.savez(fname, Q=Q, r=r, mlist=np.array(mlist))
        return Q
