from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
    Extension("thermal_analysis",  ["thermal_analysis.pyx"]),
    #   ... all your modules that need be compiled ...
                                                                ]
setup(
    name = 'Thermal Analysis',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
