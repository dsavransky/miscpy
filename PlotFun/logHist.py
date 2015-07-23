import matplotlib.pyplot as plt
import numpy as np

def logHist(X, N=30,fig=None, noclear=False, pdf=False, **kywds):
    '''
    Plot logarithmic histogram or probability density function from
    sampled data.

    Args:
        X (numpy.ndarray): 1-D array of sampled values
        N (Optional[int]): Number of bins (default 30)
        fig (Optional[int]): Figure number (default None)
        noclear (Optioanl[bool]): Clear figure (default False)
        pdf (Optional[bool]): If True normalize by bin width (default False)
            and display as curve instead of bar chart.
            Note: results are always normalized by number of samples
        **kywds: Arbitrary keyword arguments passed to matplotlib.pyplot.bar 
            (or matplotlib.pyplot.semilogx if pdf is True)

    Returns:
        x (ndarray): abscissa values of frequencies
        n (ndarray): (normalized) frequency values
    
    '''

    x = np.logspace(np.log10(np.min(X)),np.log10(np.max(X)),N+1)
    n,x = np.histogram(X,bins=x)

    n = n/float(X.size)

    plt.figure(fig)
    if not noclear: plt.clf()
    if pdf:
        n /= np.diff(x)
        x = x[:-1]+np.diff(x)/2
        plt.semilogx(x,n,**kywds)
    else:
        plt.bar(x[:len(x)-1],n,width=np.diff(x),**kywds)
    
    a = plt.gca()
    a.set_xlim(10.**np.floor(np.log10(np.min(X))),10.**np.ceil(np.log10(np.max(X))))
    a.set_xscale('log')
    plt.axis()
    
    return x,n


