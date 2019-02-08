import numpy
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension(
    name="CyNbody3",
    sources=["CyNbody3.pyx", "nbodyRKN1210_c.c"],
    include_dirs = [numpy.get_include()],
    )]

setup(
    name = 'CyNbody3',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
)
