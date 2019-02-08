from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

ext_modules = [Extension("CyRKN1210",     ["CyRKN1210.pyx"], include_dirs=[numpy.get_include()]),
               Extension("CyNbodySystem", ["CyNbodySystem.pyx"], include_dirs=[numpy.get_include()] ),
               Extension("CNbodySystem",  ["CNbodySystem.pyx", "nbodyC.c"],include_dirs = [numpy.get_include()]),
               Extension("CNbodyRKN1210", ["CNbodyRKN1210.pyx", "nbodyRKN1210_c.c"],include_dirs = [numpy.get_include()])
              ]

setup(
  name = 'Nbody RKN tools',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
