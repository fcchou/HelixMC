######
README
######

HelixMC
=======

HelixMC is a software package for Monte-Carlo (MC) simulations of DNA/RNA
helices using the base-pair level model, coded with Python. HelixMC is
distributed under the GPLv3 licence.

The project is authored by Fang-Chieh Chou in 2013, under the supervison of Dr. Rhiju Das, at the Biochemistry Department of Stanford Unviersity.

Links
=====

- Source code: https://github.com/fcchou/HelixMC
- HTML documentation: http://fcchou.github.com/HelixMC/
- Das Lab @ Stanford: http://www.stanford.edu/~rhiju/

Dependencies
============

The required dependencies to build the software are Python >= 2.7,
Numpy >= 1.6, SciPy >= 0.10, Matplotlib >= 1.1.0, and a working C/C++ compiler.

Install
=======

This package uses distutils, which is the default way of installing
python modules. To install in your home directory, use::

  python setup.py install --home

To install for all users on Unix/Linux::

  python setup.py build
  sudo python setup.py install
