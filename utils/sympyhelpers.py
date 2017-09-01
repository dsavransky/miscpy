'''
sympy helper functions
'''

import sympy
from sympy import *

th,ph,psi,thd,phd,psid,thdd,phdd,psidd = symbols('theta,phi,psi,thetadot,phidot,psidot,thetaddot,phiddot,psiddot')
w1,w2,w3 = symbols('omega_1,omega_2,omega_3')
t,g,m,h = symbols('t,g,m,h')
circmat = Matrix([eye(3)[2-j,:] for j in range(3)]) #define circulant matrix
polarframe = ['\mathbf{e_r}','\mathbf{e_\\theta}','\mathbf{e_z}']
sphericalframe = ['\mathbf{e_\\phi}','\mathbf{e_\\theta}','\mathbf{e_\\rho}']


def difftotal(expr, diffby, diffmap):
    """Take the total derivative with respect to a variable.
    Originally by Chris Wagner (http://robotfantastic.org/total-derivatives-in-sympy.html)


    Example:

        theta, t, theta_dot = symbols("theta t theta_dot")
        difftotal(cos(theta), t, {theta: theta_dot})

    returns

        -theta_dot*sin(theta)
    """
    # Replace all symbols in the diffmap by a functional form
    fnexpr = expr.subs({s:s(diffby) for s in diffmap})
    # Do the differentiation
    diffexpr = diff(fnexpr, diffby)
    # Replace the Derivatives with the variables in diffmap
    derivmap = {Derivative(v(diffby), diffby):dv 
                for v,dv in diffmap.iteritems()}
    finaldiff = diffexpr.subs(derivmap)
    # Replace the functional forms with their original form
    return finaldiff.subs({s(diffby):s for s in diffmap})

def difftotalmat(mat,diffby,diffmap):
    """Apply total derivative element by element to a matrix"""
    return Matrix([difftotal(x,diffby,diffmap) for x in mat]).reshape(*mat.shape)


def transportEq(vec, diffby, diffmap, omega):
    """Apply the transport equation to a vector:
    Frame I d/dt(vec) = Frame B d/dt(vec) + omega \cross vec
    Assumes that vector vec and angular velocity vector omega
    are given in terms of Frame B unit vectors (in the same, 
    right-handed order, i.e. b_1, b_2, b_3). diffby is the 
    variable with respect to which we differentiate (typically 
    time) and diffmap is a dictionary of derivatives of dependent
    variables. See difftotal for details."""

    return difftotalmat(vec,diffby,diffmap) + omega.cross(vec)


def rotMat(axis,angle):
    """ Returns the DCM ({}^B C^A) for a frame rotation of angle about 
    the specified axis     
    """
    if axis == 1:
        return Matrix(([1,0,0],[0,cos(angle),sin(angle)],[0,-sin(angle),cos(angle)]))
    elif axis == 2:
        return Matrix(([cos(angle),0,-sin(angle)],[0,1,0],[sin(angle), 0, cos(angle)]))
    elif axis == 3:
        return Matrix(([cos(angle),sin(angle),0],[-sin(angle),cos(angle),0],[0,0,1]))
    else:
        return -1

def parallelAxis(I_G,r_QG,m):
    """ Applies the parallel axis theorem to matrix of inertia I_G 
    (assumed to be about the COM) to find the matrix of inertia I_Q
    where the vector from G to Q is r_QG and the total mass of the 
    body is m.  I_G and I_QG are assumed to be with respect to the 
    same body frame."""
    
    return I_G + m*((r_QG.transpose()*r_QG)[0]*eye(3) - r_QG*r_QG.transpose())

def skew(v):
    """Given 3x1 vector v, return skew-symmetric matrix"""

    assert ((hasattr(v,'__iter__') or isinstance(v,Matrix)) and len(v)==3),\
            "v must be an iterable of length 3."
    
    return Matrix([
        [0,   -v[2], v[1]],
        [v[2],  0,  -v[0]],
        [-v[1],v[0],  0  ]
    ])

def calcDCM(n,th):
    """Calculated the DCM  ({}^A C^B) for a rotation of angle theta about an
    axis n"""

    return eye(3)*cos(th)+skew(n)*sin(th)+(1-cos(th))*n*n.T

def DCM2angVel(dcm,diffmap):
    """ Given a direction cosine matrix (DCM) transforming from
    frame A to frame B ({}^B C^A) returns the angular velocity of
    frame B in frame A ({}^A\omega^B). The returned vector will be 
    the components in frame B.

    Assumes that differentiation is with respect to symbol t.  
    diffmap is defined as in difftotal
    """
    
    #s = solve(dcm*difftotalmat(dcm,t,diffmap).T-skew([w1,w2,w3]),(w1,w2,w3))
    #return Matrix([s[w1],s[w2],s[w3]])

    tmp = dcm*difftotalmat(dcm,t,diffmap).T
    return simplify(Matrix([tmp[2,1],tmp[0,2],tmp[1,0]]))


def mat2vec(mat,basis='e'):
    """ Transform matrix representation of a vector to the vector equation
    for a given basis.

    basis can be either a string (i.e. 'e' becomes basis e_1,e_2,e3) or an
    iterable of length 3.

    mat is assumed to be sympy matrix representing a column vector
    """

    assert isinstance(basis,str) or (hasattr(basis,'__iter__') and len(basis)==3),\
            "basis input must be a string or iterable of length 3."

    if isinstance(basis,str):
        basis = ['\mathbf{'+basis+'}_'+str(j) for j in range(1,4)]

    basissyms = symbols(basis,commutative=False)
    basisvec = Matrix(basissyms)

    return (mat.T*basisvec)[0]

def fancyMat(prefix,shape):
    """ Create an indexed matrix akin to symarray using the prefix
    and with dimensions given by shape (this must be a 2 element iterable)

    Indexing is 1-based.

    Example:
        fancyMat('{}^\mathcal{B}C^{\mathcal{A}}',(3,3))
    """
    
    M = []
    for r in range(1,shape[0]+1):
        row = []
        for c in range(1,shape[1]+1):
            row.append(prefix+'_{'+str(r)+str(c)+'}')
        M.append(row)
    M = Matrix(symbols(M))

    return M

def EulerAngSet(rots, angs):
    """ Calculate the total DCM for an Euler Angle set defined by 
    ordered rotations about body-fixed axes rots of angles angs.
    
    Final DCM is {}^B C^A where the first rotation and angle in rots
    and angs refer to rotation starting from frame A (and the final
    entry in rots and angs is rotation into frame B).
    """

    assert (hasattr(rots,'__iter__') and len(rots)==3),\
            "rots must be an iterable of length 3."

    assert (hasattr(angs,'__iter__') and len(angs)==3),\
            "v must be an iterable of length 3."

    DCM = eye(3)
    for rot,ang in zip(rots,angs):
        DCM = rotMat(rot,ang)*DCM

    return simplify(DCM)

def EulerLagrange(L, qs, diffmap):
    """ Apply the Euler-Lagrange equations to Lagrangian L, with
    generalized coordinates qs, with their derivatives identified by the
    differentation map diffmap (see difftotal).

    Note that qs must be iterable even if this is only a one dof system.

    diffmap must map all qs to their second derivatives (via the first 
    derivative).  i.e., for any q = \theta, diffmap must contain \dot\theta
    and \ddot\theta.

    The returned value will be a system of equations for the second 
    time derivatives of the generalized coordinates.
    """

    assert hasattr(qs,'__iter__'),\
            "qs must be an iterable."

    eqs = [simplify(difftotal(L.diff(diffmap[q]),t,diffmap) - L.diff(q)) \
            for q in qs]

    return simplify(solve(eqs,[diffmap[diffmap[q]] for q in qs]))

