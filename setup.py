#! /usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
import helixmc
import numpy

DISTNAME = 'helixmc'
DESCRIPTION = "Python-based Monte Carlo simulator for DNA/RNA duplexes."
LONG_DESCRIPTION = open('README.rst').read()
AUTHOR = 'Fang-Chieh Chou'
AUTHOR_EMAIL = 'fcchou1986@gmail.com'
URL = 'http://fcchou.github.com/HelixMC/'
LICENSE = 'GPL'
DOWNLOAD_URL = 'https://github.com/fcchou/HelixMC/tarball/master'
VERSION = helixmc.__version__
PLATFORMS = ["Linux", "Mac OS-X", "Unix"]
ext_modules = [Extension(
    "helixmc._util_cython", ["helixmc/_util_cython.c"],
    include_dirs=[numpy.get_include()]
)]
package_data = {'helixmc': ['data/*.npz', 'data/*.npy']}

setup(
    name=DISTNAME,
    author=AUTHOR,
    author_email=AUTHOR_EMAIL,
    description=DESCRIPTION,
    license=LICENSE,
    url=URL,
    version=VERSION,
    platforms=PLATFORMS,
    download_url=DOWNLOAD_URL,
    long_description=LONG_DESCRIPTION,
    packages=['helixmc', 'helixmc.tests'],
    scripts=['helixmc-run'],
    package_data=package_data,
    ext_modules=ext_modules,
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
