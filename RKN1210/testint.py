from scipy.io import loadmat, savemat
import numpy as np

#load data
data = loadmat('solarSystemData.mat')
mus = np.squeeze(data['mus'])
p0 = np.squeeze(data['p0'])
v0 = np.squeeze(data['v0'])

import CNbodyRKN1210, CNbodySystem, CyNbodySystem, CyRKN1210, NbodySystem, RKN1210

py_nbs = NbodySystem.NbodySystem(mus)
cy_nbs = CyNbodySystem.NbodySystem(mus)
c_nbs = CNbodySystem.NbodySystem(mus)


tin = np.linspace(0,10*365.25,1000)
#python integrator using python nbody eq
%time t1,y1,dy1 = RKN1210.RKN1210(py_nbs.nbody_eq,tin,p0,v0)

#python integrator using python vectorized nbody eq
%time t2,y2,dy2 = RKN1210.RKN1210(py_nbs.nbody_eq_vec,tin,p0,v0)

#python integrator using cythonized nbody_eq 
%time t3,y3,dy3 = RKN1210.RKN1210(cy_nbs.nbody_eq,tin,p0,v0)

#python integrator using cythonized C nbody_eq 
%time t4,y4,dy4 = RKN1210.RKN1210(c_nbs.nbody_eq,tin,p0,v0)

#cython integrator using python nbody eq
%time t5,y5,dy5 = CyRKN1210.RKN1210(py_nbs.nbody_eq,tin,p0,v0)

#cython integrator using python vectorized nbody eq
%time t6,y6,dy6 = CyRKN1210.RKN1210(py_nbs.nbody_eq_vec,tin,p0,v0)

#cython integrator using cythonized nbody_eq 
%time t7,y7,dy7 = CyRKN1210.RKN1210(cy_nbs.nbody_eq,tin,p0,v0)

#cython integrator using cythonized C nbody_eq 
%time t8,y8,dy8 = CyRKN1210.RKN1210(c_nbs.nbody_eq,tin,p0,v0)


#everything in pure C with cython wrapping
%time t9,y9,dy9 = CNbodyRKN1210.nbodyRKN1210(tin, p0, v0, mus)
