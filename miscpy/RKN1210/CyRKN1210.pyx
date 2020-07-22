import sys
import numpy as np
cimport numpy as np
DTYPE = np.double
ctypedef np.double_t DTYPE_t
cimport cython
@cython.boundscheck(False)
@cython.wraparound(False)

def RKN1210(fun, tspan0, np.ndarray[DTYPE_t, ndim=1] y0, np.ndarray[DTYPE_t, ndim=1] dy0, double abstol=1e-14):
    """RKN1210
    
    12th/10th order Runge—Kutta—Nyström method for the numerical integration of
    second-order differential equations of the form u''(t) = f(t, u)

    Args:         
        fun (callable):
            Function to be integrated.  Must take inputs t,y for 
            scalar time t and state vector y and return y''.
        tspan (scalar or array):
            Times to find solutions at.  If scalar, assumes integration
            is between 0 and tspan.
        y0 (array)L
            Initial state.
        dy0 (array):
            Initial state derivative (y').  Must be of same dimension as y0
        abstol (scalar float):
            Optional absolute tolerance for integration convergence.  Defaults
            to 1e-14.  This controls the adaptive step-size of the internal 
            time step used.
    
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

        Equivalent to RKN1210 but fully cythonized!
    
    """

    
    cdef np.ndarray[DTYPE_t, ndim=1] c = np.array([0.0e0,
        2.0e-2,
        4.0e-2,
        1.0e-1,
        (1.+1./3)/10,
        1.6e-1,
        5.0e-2,
        2.0e-1,
        2.5e-1,
        (3.+1./3)/10.,
        5.0e-1,
        (5.+5./9)/10.,
        7.5e-1,
        8.57142857142857142857142857143e-1,
        9.45216222272014340129957427739e-1,
        1.0e0,
        1.0e0],dtype=DTYPE)

    '''
    matrix A is lower triangular. It's easiest to 
    load the coefficients row-by-row:
    '''   
    cdef np.ndarray[DTYPE_t, ndim=2] A = np.zeros((17,17), dtype=DTYPE)
    A[1,0] = 2.0e-4  

    A[2,0:2] = np.array([2.66666666666666666666666666667e-4,   
                5.33333333333333333333333333333e-4],dtype=DTYPE)
            
    A[3,0:3] = np.array([2.91666666666666666666666666667e-3,
                -4.16666666666666666666666666667e-3,
                6.25e-3],dtype=DTYPE)    
            
    A[4,0:4] = np.array([1.64609053497942386831275720165e-3,
                0.0e0,
                5.48696844993141289437585733882e-3,
                1.75582990397805212620027434842e-3],dtype=DTYPE)   
            
    A[5,0:5] = np.array([1.9456e-3,
                0.0e0,
                7.15174603174603174603174603175e-3,
                2.91271111111111111111111111111e-3,
                7.89942857142857142857142857143e-4],dtype=DTYPE)    
            
    A[6,0:6] = np.array([5.6640625e-4,
                0.0e0,
                8.80973048941798941798941798942e-4,
                -4.36921296296296296296296296296e-4,
                3.39006696428571428571428571429e-4,
                -9.94646990740740740740740740741e-5],dtype=DTYPE)    
            
    A[7,0:7] = np.array([3.08333333333333333333333333333e-3,
                0.0e0,
                0.0e0,
                1.77777777777777777777777777778e-3,
                2.7e-3,
                1.57828282828282828282828282828e-3,
                1.08606060606060606060606060606e-2],dtype=DTYPE)    
            
    A[8,0:8] = np.array([3.65183937480112971375119150338e-3,
                0.0e0,
                3.96517171407234306617557289807e-3,
                3.19725826293062822350093426091e-3,
                8.22146730685543536968701883401e-3,
                -1.31309269595723798362013884863e-3,
                9.77158696806486781562609494147e-3,
                3.75576906923283379487932641079e-3],dtype=DTYPE)

    A[9,0:9] = np.array([3.70724106871850081019565530521e-3,
                0.0e0,
                5.08204585455528598076108163479e-3,
                1.17470800217541204473569104943e-3,
                -2.11476299151269914996229766362e-2,
                6.01046369810788081222573525136e-2,
                2.01057347685061881846748708777e-2,
                -2.83507501229335808430366774368e-2,
                1.48795689185819327555905582479e-2],dtype=DTYPE)

    A[10,0:10] = np.array([3.51253765607334415311308293052e-2,
                0.0e0,
                -8.61574919513847910340576078545e-3,
                -5.79144805100791652167632252471e-3,
                1.94555482378261584239438810411e0,
                -3.43512386745651359636787167574e0,
                -1.09307011074752217583892572001e-1,
                2.3496383118995166394320161088e0,
                -7.56009408687022978027190729778e-1,
                1.09528972221569264246502018618e-1],dtype=DTYPE)

    A[11,0:11] = np.array([2.05277925374824966509720571672e-2,
                0.0e0,
                -7.28644676448017991778247943149e-3,
                -2.11535560796184024069259562549e-3,
                9.27580796872352224256768033235e-1,
                -1.65228248442573667907302673325e0,
                -2.10795630056865698191914366913e-2,
                1.20653643262078715447708832536e0,
                -4.13714477001066141324662463645e-1,
                9.07987398280965375956795739516e-2,
                5.35555260053398504916870658215e-3],dtype=DTYPE)

    A[12,0:12] = np.array([-1.43240788755455150458921091632e-1,
                0.0e0,
                1.25287037730918172778464480231e-2,
                6.82601916396982712868112411737e-3,
                -4.79955539557438726550216254291e0,
                5.69862504395194143379169794156e0,
                7.55343036952364522249444028716e-1,
                -1.27554878582810837175400796542e-1,
                -1.96059260511173843289133255423e0,
                9.18560905663526240976234285341e-1,
                -2.38800855052844310534827013402e-1,
                1.59110813572342155138740170963e-1],dtype=DTYPE)

    A[13,0:13] = np.array([8.04501920552048948697230778134e-1,
                0.0e0,
                -1.66585270670112451778516268261e-2,
                -2.1415834042629734811731437191e-2,
                1.68272359289624658702009353564e1,
                -1.11728353571760979267882984241e1,
                -3.37715929722632374148856475521e0,
                -1.52433266553608456461817682939e1,
                1.71798357382154165620247684026e1,
                -5.43771923982399464535413738556e0,
                1.38786716183646557551256778839e0,
                -5.92582773265281165347677029181e-1,
                2.96038731712973527961592794552e-2],dtype=DTYPE)

    A[14,0:14] = np.array([-9.13296766697358082096250482648e-1,
                0.0e0,
                2.41127257578051783924489946102e-3,
                1.76581226938617419820698839226e-2,
                -1.48516497797203838246128557088e1,
                2.15897086700457560030782161561e0,
                3.99791558311787990115282754337e0,
                2.84341518002322318984542514988e1,
                -2.52593643549415984378843352235e1,
                7.7338785423622373655340014114e0,
                -1.8913028948478674610382580129e0,
                1.00148450702247178036685959248e0,
                4.64119959910905190510518247052e-3,
                1.12187550221489570339750499063e-2],dtype=DTYPE)

    A[15,0:15] = np.array([-2.75196297205593938206065227039e-1,
                0.0e0,
                3.66118887791549201342293285553e-2,
                9.7895196882315626246509967162e-3,
                -1.2293062345886210304214726509e1,
                1.42072264539379026942929665966e1,
                1.58664769067895368322481964272e0,
                2.45777353275959454390324346975e0,
                -8.93519369440327190552259086374e0,
                4.37367273161340694839327077512e0,
                -1.83471817654494916304344410264e0,
                1.15920852890614912078083198373e0,
                -1.72902531653839221518003422953e-2,
                1.93259779044607666727649875324e-2,
                5.20444293755499311184926401526e-3],dtype=DTYPE)

    A[16,0:16] = np.array([1.30763918474040575879994562983e0,
                0.0e0,
                1.73641091897458418670879991296e-2,
                -1.8544456454265795024362115588e-2,
                1.48115220328677268968478356223e1,
                9.38317630848247090787922177126e0,
                -5.2284261999445422541474024553e0,
                -4.89512805258476508040093482743e1,
                3.82970960343379225625583875836e1,
                -1.05873813369759797091619037505e1,
                2.43323043762262763585119618787e0,
                -1.04534060425754442848652456513e0,
                7.17732095086725945198184857508e-2,
                2.16221097080827826905505320027e-3,
                7.00959575960251423699282781988e-3,
                0.0e0],dtype=DTYPE)  
      
    A = A.T           

    #Bhat (high-order b)
    cdef np.ndarray[DTYPE_t, ndim=1] Bhat = np.array([1.21278685171854149768890395495e-2,
            0.0e0,
            0.0e0,
            0.0e0,
            0.0e0,
            0.0e0,
            8.62974625156887444363792274411e-2,
            2.52546958118714719432343449316e-1,
            -1.97418679932682303358307954886e-1,
            2.03186919078972590809261561009e-1,
            -2.07758080777149166121933554691e-2,
            1.09678048745020136250111237823e-1,
            3.80651325264665057344878719105e-2,
            1.16340688043242296440927709215e-2,
            4.65802970402487868693615238455e-3,
            0.0e0,
            0.0e0], dtype=DTYPE)

    #BprimeHat (high-order b-prime)
    cdef np.ndarray[DTYPE_t, ndim=1] Bphat = np.array([1.21278685171854149768890395495e-2,
            0.0e0,
            0.0e0,
            0.0e0,
            0.0e0,
            0.0e0,
            9.08394342270407836172412920433e-2,
            3.15683697648393399290429311645e-1,
            -2.63224906576909737811077273181e-1,
            3.04780378618458886213892341513e-1,
            -4.15516161554298332243867109382e-2,
            2.46775609676295306562750285101e-1,
            1.52260530105866022937951487642e-1,
            8.14384816302696075086493964505e-2,
            8.50257119389081128008018326881e-2,
            -9.15518963007796287314100251351e-3,
            2.5e-2], dtype=DTYPE)

    #B (low-order b)
    cdef np.ndarray[DTYPE_t, ndim=1] B = np.array([1.70087019070069917527544646189e-2,
            0.0e0,
            0.0e0,
            0.0e0,
            0.0e0,
            0.0e0,
            7.22593359308314069488600038463e-2,
            3.72026177326753045388210502067e-1,
            -4.01821145009303521439340233863e-1,
            3.35455068301351666696584034896e-1,
            -1.31306501075331808430281840783e-1,
            1.89431906616048652722659836455e-1,
            2.68408020400290479053691655806e-2,
            1.63056656059179238935180933102e-2,
            3.79998835669659456166597387323e-3,
            0.0e0,
            0.0e0], dtype=DTYPE)

    #Bprime (low-order bprime)
    cdef np.ndarray[DTYPE_t, ndim=1] Bp = np.array([1.70087019070069917527544646189e-2,
            0.0e0,
            0.0e0,
            0.0e0,
            0.0e0,
            0.0e0,
            7.60624588745593757356421093119e-2,
            4.65032721658441306735263127583e-1,
            -5.35761526679071361919120311817e-1,
            5.03182602452027500044876052344e-1,
            -2.62613002150663616860563681567e-1,
            4.26221789886109468625984632024e-1,
            1.07363208160116191621476662322e-1,
            1.14139659241425467254626653171e-1,
            6.93633866500486770090602920091e-2,
            2.0e-2,
            0.0e0], dtype=DTYPE)

    #format time input
    cdef np.ndarray[DTYPE_t, ndim = 1] tspan
    if isinstance(tspan0, (list,np.ndarray)):
        tspan = np.array(tspan0,dtype=DTYPE)
        userTime = True
    elif isinstance(tspan0, (int,long,float)):
        tspan = np.array([0.0,tspan0],dtype=DTYPE)
        userTime = False
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
    if (y0.size != dy0.size):
        raise Exception("y0 and dy0 must be n x 1 (or 1 x n) vectors of the same dimension.")   
    assert y0.dtype == DTYPE and dy0.dtype == DTYPE
        
    #initialize
    cdef np.ndarray[DTYPE_t, ndim=1] y = y0
    cdef np.ndarray[DTYPE_t, ndim=1] dy = dy0
    cdef double pow = 1./12
    cdef double t = t0
    cdef np.ndarray[DTYPE_t, ndim=2] f = np.zeros((y.size,17), dtype=DTYPE)
    cdef double hmin = abs(tf - t0)/1e12
    cdef double hmax = tf - t0
    try:
        f[:,0] = fun(t,y)
    except:
        e = sys.exc_info()[1]
        print "Error %s" % e

    #default initial step
    cdef double h = (abstol ** pow) / max(max(abs(np.hstack((dy, f[:, 0])))), 1e-4)
    h = min(hmax,max(h,hmin))
    h = tdir*h
    cdef double h0

    #pre-allocate output arrays
    cdef np.ndarray[DTYPE_t, ndim=1] thist
    cdef np.ndarray[DTYPE_t, ndim=2] yhist, dyhist
    cdef int histSize
    if userTime:
        thist = np.zeros(tspan.size, dtype=DTYPE)
        yhist = np.zeros((y.size,tspan.size), dtype=DTYPE)
        dyhist = np.zeros((dy.size,tspan.size), dtype=DTYPE)
    else:
        histSize = max(round(abs(tf/h)),2)
        thist = np.zeros(histSize, dtype=DTYPE)
        yhist = np.zeros((y.size,histSize), dtype=DTYPE)
        dyhist = np.zeros((y.size,histSize), dtype=DTYPE)

    #write ICs to history
    thist[0] = t0
    yhist[:,0] = y
    dyhist[:,0] = dy

    #pre-allocate things used in loop
    cdef int counter = 1
    cdef double tmp, delta
    cdef np.ndarray[DTYPE_t, ndim=1] fBphat = np.zeros(y.size , dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] fBhat = np.zeros(y.size , dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] fB = np.zeros(y.size , dtype=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] fBp = np.zeros(y.size , dtype=DTYPE)
    
    #enter main loop
    while (abs(t - tf) > 0):
        h0 = h
        #check to see if this step will be written to history
        if not userTime:
            addToHist = True
            #expand arrays as necessary
            if counter == thist.size:
                tmp = round(abs((tf-t)/np.median(thist[1:] - thist[:thist.size-1])))
                yhist = np.hstack((yhist,np.zeros((y.size,tmp), dtype=DTYPE)))
                dyhist = np.hstack((dyhist,np.zeros((dy.size,tmp), dtype=DTYPE)))
                thist = np.hstack((thist,np.zeros(tmp, dtype=DTYPE)))

            #if going past final time, fix time step
            if (tdir*(t + h - tf) > 0): h = tf - t
            
        else:
            addToHist = False
            #if this step goes past a usertime decrease it and add to hist
            if (tdir*(t + h - tspan[counter]) >= -hmin):
                addToHist = True
                h = tspan[counter] - t

        #compute 2nd derivative
        for j in range(17):
            f[:,j] = fun( t + c[j]*h, y + c[j]*h*dy + (h ** 2.)*np.dot(f,A[:,j]) )

        if np.isnan(f).any() or np.isinf(f).any():
            raise Warning("Second derivative computed with infinite or NaN values")
            #return thist[:counter],yhist[:,:counter],dyhist[:,:counter]

        #estimate the error
        #fBphat = np.dot(f,Bphat)
        #fBhat  = np.dot(f,Bhat)
        #fB = np.dot(f,B)
        #fBp = np.dot(f,Bp)
        
        dotAb(f,Bphat,fBphat)
        dotAb(f,Bhat,fBhat)
        dotAb(f,B,fB)
        dotAb(f,Bp,fBp)
        delta = max( max(abs( (h ** 2.)*(fBhat - fB))), max(abs(h*(fBphat - fBp))) ) * h

        #if error is ok, update the solution & write to history if needed
        if ( delta <= abstol ):
                t = t + h
                y = y + h*dy + ( h ** 2.) *fBhat
                dy = dy + h*fBphat

                if addToHist:
                    thist[counter] = t
                    yhist[:,counter] = y
                    dyhist[:,counter] = dy
                    counter += 1
        
        #update next step size
        if (delta != 0):
            h = tdir * min(hmax, 0.9*abs(h0)*((abstol/delta) ** pow))
            if (abs(h) < hmin):
                raise Warning("Step size below minimum at time "+repr(t)+". Singularity likely.")
                #return thist[:counter],yhist[:,:counter],dyhist[:,:counter]

    #finished loop so return
    return thist[:counter],yhist[:,:counter],dyhist[:,:counter]
        

def dotAb(np.ndarray[DTYPE_t, ndim=2] A, np.ndarray[DTYPE_t, ndim=1] b, np.ndarray[DTYPE_t, ndim=1] out = None):

    cdef Py_ssize_t i, j, k
    cdef DTYPE_t s
    #out = np.zeros(A.shape[0], dtype=DTYPE)
    for i in range(A.shape[0]):
        s = 0
        for k in range(A.shape[1]):
            s +=A[i,k]*b[k]
        out[i] = s
