'''
sympy helper functions
'''

from sympy import *

th,ph,psi,thd,phd,psid,thdd,phdd,psidd = symbols('theta,phi,psi,thetadot,phidot,psidot,thetaddot,phiddot,psiddot')
w1,w2,w3 = symbols('omega_1,omega_2,omega_3')
t,g,m,h = symbols('t,g,m,h')

#By Chris Wagner (http://robotfantastic.org/total-derivatives-in-sympy.html):
def difftotal(expr, diffby, diffmap):
    """Take the total derivative with respect to a variable.

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
    """ Returns the DCM for a frame rotation of angle about 
    the specified axis ({}^B C^A)
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
    
    return I_G + m*((r_QG.transpose()*r_QG)[0]*eye(3) + r_QG*r_QG.transpose())

def skew(v):
    """Given 3x1 vector v, return skew-symmetric matrix"""
    
    return Matrix([
        [0,   -v[2], v[1]],
        [v[2],  0,  -v[0]],
        [-v[1],v[0],  0  ]
    ])

