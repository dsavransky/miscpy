import numpy as np
from scipy.optimize import fmin_l_bfgs_b

def simpSample(f, numTest, xMin, xMax, verb=False, maxIter=1000):
    '''
    Simple Rejection Sampling

    Args:
        f (callable): Probability density function (may be 
            unnormalized).  Function or lambda function.
        numTest (int): Number of samples to generate
        xMin (float): Minimum sampling range
        xMax (float): Maximum sampling range
        verb (Optional[bool]): Print number of iterations taken
            (default False)
        maxIter (Optional[int]): Maximum number of iterations to
            take (default 1000).

    Returns:
        X (ndarray): Array of samples from given PDF

    Raises:
        ArithmeticError: if maximum number of iterations is reached
            before enough samples are collected
   
    Example:
        #sample from a log-uniform distribution
        X = simpSample(lambda x: 1./x, 1000, 0.1,100)
    '''

    #find max value
    #first do a coarse grid to get ic
    dx = np.linspace(xMin,xMax,100)
    ic = np.argmax(f(dx))
    g = lambda x: -f(x)
    M = fmin_l_bfgs_b(g,[dx[ic]],approx_grad=True,bounds=[(xMin,xMax)])
    M = f(M[0])

    #initialize
    n = 0
    X = np.zeros(numTest);
    numIter = 0;

    nSamp = np.max([2*numTest,1e6])
    while n < numTest and numIter < maxIter:
        xd = np.random.random(nSamp) * (xMax - xMin) + xMin
        yd = np.random.random(nSamp) * M
        pd = f(xd)

        xd = xd[yd < pd]
        X[n:min(n+len(xd), numTest)] = xd[:min(len(xd),numTest-n)]
        n += len(xd)
        numIter += 1

    if numIter == maxIter:
        raise ArithmeticError("Failed to converge.")

    if verb:
        print 'Finished in '+repr(numIter)+' iterations.'

    return X

