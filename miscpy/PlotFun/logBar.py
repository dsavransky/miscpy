import matplotlib.pyplot as plt
import numpy as np

def logBar(x,n,fig=None,nobar = False, noclear = False, norm=True, **kywds):

    '''
    Plot logarithmic bar plot such that all bars are evenly spaced and of 
    uniform width in log space

    Args:
        x (ndarray): left edge values of frequencies (bar left edge)
        n (ndarray): (normalized) frequency values
        fig (Optional[int]): Figure number (default None)
        noclear (Optioanl[bool]): Clear figure (default False)
        norm (Optional([bool]): Normalize frequencies by number of samples 
            (default True)
        **kywds: Arbitrary keyword arguments passed to matplotlib.pyplot.bar 
            (or matplotlib.pyplot.semilogx if pdf is True)

    Returns:
        x (ndarray): abscissa values of frequencies
        n (ndarray): (normalized) frequency values
    
    '''

    
    plt.figure(fig)
    if not noclear: plt.clf()
    if nobar:
        plt.semilogx(x[1:],n/float(np.sum(n)),**kywds)
    else:
        plt.bar(x[:len(x)-1],n/float(np.sum(n)),width=np.diff(x),**kywds)

    a = plt.gca()
    a.set_xlim(10.**np.floor(np.log10(min(x))),10.**np.ceil(np.log10(max(x))))
    a.set_xscale('log')
    plt.axis()


