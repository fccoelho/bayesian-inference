import cython
cimport numpy as np

cdef class Sampler:
    @cython.locals(D=np.ndarray, Dbar=np.ndarray, i=cython.int)
    cpdef double DIC():
    cpdef _accept(self, double last_lik, double lik)
    @cython.locals(i=cython.int, c=cython.int, step=cython.int)
    cpdef _propose(self,parpriors=[],  parnames = [])
