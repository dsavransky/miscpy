import numpy as np

class NbodySystem:
    """NbodySystem
    
    This class generates an object for calculating the mutual accelerations due to
    Newtonian gravity in an N-body system.  The object provides two methods for 
    performing the calculation: one (nbody_eq) based on a standard looping algorithm, 
    and a second (nobdy_eq_vec), vectorized variant.

    Args:         
        mus (iterable):
            user specified values for the gravitational parameter of each body

    
    """

    def __init__(self, mus):
        n = mus.size
        self.inds = np.reshape(np.arange(3*n),(n,3)).T
        cols = np.squeeze(np.tile(np.arange(n),(1,n)))
        self.ks = cols[np.nonzero(np.reshape(np.ones((n,n)) - np.eye(n),n ** 2))]
        self.js = np.reshape(np.tile(np.arange(n),(n-1,1)).T,n*(n-1))
        self.inds2 = np.reshape(np.reshape(np.arange(n*(n-1)),(n,n-1)).T,n*(n-1))
        self.mus = np.squeeze(mus)
        self.n = n

    def nbody_eq_vec(self,t,x):
        '''
        Vectorized N-body equations.

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

        n = self.n
        
        #r_k/j = r_k/o - r_j/o
        rkj = x[self.inds[:,self.ks]] - x[self.inds[:,self.js]]

        #d^2/dt^2(r_k/j) = G*sum_{k ~= j}(mu_k*r_k/j)/(r_k/j^3)
        rkj = np.dot(rkj,np.dot(np.diag(1/np.sqrt(sum(rkj ** 2.,0)) ** 3),np.diag(self.mus[self.ks])))

        return sum(np.reshape(rkj[:,self.inds2].T,(n-1,3*n)),0)
        

    def nbody_eq2(self, t, x):
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

        dx = np.zeros(x.size)
        
        for j in range(self.n):
            temp = np.zeros(3)
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
