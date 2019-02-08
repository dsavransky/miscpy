import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension(
    name="CyNbody2",
    sources=["CyNbody2.pyx", "nbodyC.c"],
    include_dirs = [numpy.get_include()],
    )]

setup(
    name = 'CyNbody2',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
)
