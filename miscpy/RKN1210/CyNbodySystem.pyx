import numpy as np
cimport numpy as np
DTYPE = np.double
ctypedef np.double_t DTYPE_t

cdef extern from "math.h":
    DTYPE_t sqrt(DTYPE_t)
cdef extern from "math.h":
    DTYPE_t pow(DTYPE_t, DTYPE_t)

cimport cython
@cython.boundscheck(False)
@cython.wraparound(False)

cdef class NbodySystem:
    """CyNbodySystem
    
    This class generates an object for calculating the mutual accelerations due to
    Newtonian gravity in an N-body system.  The object provides two methods for 
    performing the calculation: one (nbody_eq) based on a standard looping algorithm, 
    and a second (nobdy_eq_vec), vectorized variant. Identical to NbodySystem but 
    Cythonized.

    Args:         
        mus (iterable):
            user specified values for the gravitational parameter of each body

    
    """

    cdef np.ndarray mus
    cdef int n

    def __init__(self, np.ndarray[DTYPE_t, ndim=1] mus):
        self.n = mus.size
        self.mus = mus

    cdef np.ndarray[DTYPE_t, ndim = 1] _nbody_eq(self, DTYPE_t t, np.ndarray[DTYPE_t, ndim=1] x):

        cdef np.ndarray[DTYPE_t, ndim = 1] temp = np.zeros(3,dtype=DTYPE)
        cdef DTYPE_t rnorm3, mu
        cdef int j, k, i
        cdef np.ndarray[DTYPE_t, ndim = 1] dx = np.zeros(3*self.n,dtype=DTYPE)

        for j in range(self.n):
            temp = np.zeros(3,dtype=DTYPE)
            for k in range(self.n):
                mu = self.mus[k]
                if (k == j): continue
                rnorm3 = pow(pow(x[k*3] - x[j*3],2.) + pow(x[k*3+1] - x[j*3+1],2.) + pow(x[k*3+2] - x[j*3+2],2.),3./2.)
                if (rnorm3 == 0.0):
                    raise Exception("Divide by zero.")
                for i in range(3):
                    temp[i] += (x[k*3+i] - x[j*3+i])/rnorm3*mu

            for i in range(3):
                dx[j*3+i] = temp[i]

        return dx

    def nbody_eq(self, DTYPE_t t, np.ndarray[DTYPE_t, ndim=1] x):
        '''
        N-body equations.

        Calculate the mutual accelerations of a collection of N bodies defined by 
        gravitational parameters (G*M) stored in self.mu
    
        Args:
            t (scalar):
                time - dummy var required by many integrators. Not used in calculation.
            x (ndarray):
                3n x 1 vector of stacked positions or n planets where n = self.n.  
            
        Return:
            ddx (ndarray):
                Second derivative of input.

        Notes:
            Position units must be complementary (i.e., if position is AU mu must be in 
            AU^3/TU^2).

        '''

        return self._nbody_eq(t, x)
    
