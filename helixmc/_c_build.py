from distutils.core import setup
from distutils.extension import Extension
import numpy

ext_modules = [Extension(
    "_util_cython", ["_util_cython.c"], include_dirs=[numpy.get_include()])]


setup(ext_modules=ext_modules)
