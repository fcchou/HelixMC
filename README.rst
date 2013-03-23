.. -*- mode: rst -*-

HelixMC
=======

HelixMC is a software package for Monte-Carlo (MC) simulations of DNA/RNA
helices using the base-pair level model, coded with Python. HelixMC is
distributed under the GPLv3 licence.

The project is authored by Fang-Chieh Chou in 2013, under the supervison of his
research advisor Dr. Rhiju Das, at the Biochemistry Department of Stanford
Unviersity.

Links
=====

- Source code in GitHub: https://github.com/fcchou/HelixMC
- HTML documentation:

Dependencies
============

The required dependencies to build the software are Python >= 2.7,
setuptools, Numpy >= 1.6, SciPy >= 0.10 and a working C/C++ compiler.

Install
=======

This package uses distutils, which is the default way of installing
python modules. To install in your home directory, use::

  python setup.py install --home

To install for all users on Unix/Linux::

  python setup.py build
  sudo python setup.py install
