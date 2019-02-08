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
        return self._nbody_eq(t, x)
    
