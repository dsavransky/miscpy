import numpy as np
cimport numpy as np
DTYPE = np.double
ctypedef np.double_t DTYPE_t

cdef extern from "nbodyC.h":
    int nbodyC(double* dx, double* x, double* mus, int n )

cdef class NbodySystem:
    """CNbodySystem
    
    This class generates an object for calculating the mutual accelerations due to
    Newtonian gravity in an N-body system.  The object provides two methods for 
    performing the calculation: one (nbody_eq) based on a standard looping algorithm, 
    and a second (nobdy_eq_vec), vectorized variant. Identical to NbodySystem but 
    Cythonized wrapper around pure C implementation.

    Args:         
        mus (iterable):
            user specified values for the gravitational parameter of each body

    
    """

    cdef np.ndarray mus
    cdef int n

    def __init__(self, np.ndarray[DTYPE_t, ndim=1] mus):
        self.n = mus.size
        self.mus = mus

    def nbody_eq(self, double t, np.ndarray[DTYPE_t,ndim=1] x):
        cdef np.ndarray[DTYPE_t, ndim=1] dx = np.zeros(x.size,dtype=DTYPE)
        nbodyC( <double*> dx.data, <double*> x.data, <double*> self.mus.data, self.n )
        return dx
