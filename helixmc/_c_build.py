from distutils.core import setup
from distutils.extension import Extension

ext_modules = [Extension("_util_cython", ["_util_cython.c"])]

setup(
  ext_modules = ext_modules
)
