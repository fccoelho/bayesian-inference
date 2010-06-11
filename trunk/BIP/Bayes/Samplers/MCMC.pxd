import cython
cimport numpy as np

cdef class Sampler:
    @cython.locals(D=np.ndarray, Dbar=np.ndarray, i=cython.int)
    cpdef double DIC():
    cpdef _accept(self, double last_lik, double lik)
    @cython.locals(i=cython.int, c=cython.int, step=cython.int)
    cpdef _propose(self,parpriors=[],  parnames = [])

cdef class Dream:
    @cython.locals(alpha=cython.double)
    cpdef double _alpha1(self, double p1, double p2)
    
    @cython.locals(CR=cython.double, delta=cython.int,gam=cython.double, d=cython.int)
    cpdef _chain_evolution(self, proptheta, propphi, pps, liks)
    
    @cython.locals(means=np.ndarray, q1=cython.double, q3=cython.double, iqr=cython.double)
    cpdef _det_outlier_chains(self, int step)
    
    @cython.locals()
    cpdef _get_post_prob(self, theta, prop)
