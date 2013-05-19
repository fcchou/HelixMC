#! /usr/bin/env python

import os
from distutils.core import setup
from distutils.extension import Extension
import helixmc
import numpy

DISTNAME = 'helixmc'
DESCRIPTION = "Python-based Monte Carlo simulator for DNA/RNA helices."
LONG_DESCRIPTION = open('README.rst').read()
AUTHOR = 'Fang-Chieh Chou'
AUTHOR_EMAIL = 'fcchou@stanford.edu'
URL = 'http://fcchou.github.com/HelixMC/'
LICENSE = 'GPL'
DOWNLOAD_URL = 'https://github.com/fcchou/HelixMC'
VERSION = str(helixmc.__version__)
ext_modules = [ Extension("helixmc._util_cython", ["helixmc/_util_cython.c"], include_dirs=[numpy.get_include()] ]
package_data={'helixmc': ['database/*.npz']}

setup(name=DISTNAME,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      description=DESCRIPTION,
      license=LICENSE,
      url=URL,
      version=VERSION,
      download_url=DOWNLOAD_URL,
      long_description=LONG_DESCRIPTION,
      packages=["helixmc"],
      scripts=['helixmc-run'],
      package_data=package_data,
      ext_modules = ext_modules,
      classifiers=[
          'Intended Audience :: Science/Research',
          'License :: OSI Approved',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering',
          'Operating System :: POSIX',
          'Operating System :: Unix',
          'Operating System :: MacOS'
         ],
    )
