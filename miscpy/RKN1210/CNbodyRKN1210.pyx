import numpy as np
cimport numpy as np
DTYPE = np.double
ctypedef np.double_t DTYPE_t

cdef extern from "nbodyRKN1210_c.h":
    int nbodyRKN1210_c(double* t_in, double* y_in, double* dy_in,
                       double* mus, int N, double* y, double* dy)

cimport cython
@cython.boundscheck(False)
@cython.wraparound(False)
    
def nbodyRKN1210(tspan0, np.ndarray[DTYPE_t, ndim=1] y0,
                 np.ndarray[DTYPE_t, ndim=1] dy0, np.ndarray[DTYPE_t, ndim=1] mus):
    """nbodyRKN1210
    
    12th/10th order Runge—Kutta—Nyström method for the numerical integration of
    Newtonian N-body gravitational system.

    Args:         
        tspan0 (scalar or array):
            Times to find solutions at.  If scalar, assumes integration
            is between 0 and tspan.
        y0 (array)L
            Initial state.
        dy0 (array):
            Initial state derivative (y').  Must be of same dimension as y0
        mus (array):
            Gravitational parameters 
    
    Returns:
        thist (array):
            Times corresponding to solutions entries
        yhist (2d array):
            Integrator states at each time
        dyhist (2d array):
            Integrator state dervatives at each time

    Notes:
        The construction of RKN12(10) is described in:
        High-Order Embedded Runge-Kutta-Nystrom Formulae
        J. R. DORMAND, M. E. A. EL-MIKKAWY, AND P. J. PRINCE
        IMA Journal of Numerical Analysis (1987) 7, 423-430
        
        Coefficients obtained from http://www.tampa.phys.ucl.ac.uk/rmat/test/rknint.f

        This code largely based on MATLAB implementation by Rody P.S. Oldenhuis, as
        found on Mathworks File Exchange.

        Purely a cythonized wrapper of the underlying nbodyRKN1210_c.c code.

    """
 
    #format time input
    cdef np.ndarray[DTYPE_t, ndim = 1] tspan
    if isinstance(tspan0, (list,np.ndarray)):
        tspan = np.array(tspan0,dtype=DTYPE)
    elif isinstance(tspan0, (int,long,float)):
        tspan = np.array([0.0,tspan0],dtype=DTYPE)
    else:
        raise Exception("tspan must be scallar or array.")

    cdef double t0 = tspan[0]
    cdef double tf = tspan[tspan.size-1]
    cdef double tdir = (tf - t0)/abs(tf - t0)
    if t0 == tf:
        raise Exception("t0 cannot equal tf.")
    if any((tspan[1:] - tspan[:tspan.size-1])*tdir <= 0):
        raise Exception("tspan must be monotonically increasing or decreasing.")

    #validate other inputs
    cdef int N = mus.size
    if (y0.size != dy0.size) or (y0.size != N*3):
        raise Exception("y0 and dy0 must have 3n elements and mus must have n elements.")   
    assert y0.dtype == DTYPE and dy0.dtype == DTYPE

    #initialize outputs
    cdef np.ndarray[DTYPE_t, ndim=2] yhist = np.zeros((y0.size,tspan.size), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=2] dyhist = np.zeros((dy0.size,tspan.size), dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] y_in = y0
    cdef np.ndarray[DTYPE_t, ndim=1] dy_in = dy0
    cdef np.ndarray[DTYPE_t, ndim=1] y_out = np.zeros(y0.size,dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] dy_out = np.zeros(y0.size,dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] ts

    #write ICs to history
    yhist[:,0] = y0
    dyhist[:,0] = dy0

    cdef int j, res
    for j from 1 <= j < tspan.size:
        ts = np.array((tspan[j-1],tspan[j]),dtype=DTYPE)
        res = nbodyRKN1210_c(<double*> ts.data, <double*> y_in.data, <double*>  dy_in.data,
                         <double*>  mus.data, N, <double*>  y_out.data, <double*> dy_out.data)
        if (res != 0):
            raise Exception("Integration failed.")

        y_in = y_out
        dy_in = dy_out
        yhist[:,j] = y_out
        dyhist[:,j] = dy_out
        
    return tspan, yhist, dyhist
