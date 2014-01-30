.. image:: https://travis-ci.org/fcchou/HelixMC.png?branch=master
      :target: https://travis-ci.org/fcchou/HelixMC/

######
README
######

HelixMC
=======
HelixMC is a software package for Monte-Carlo (MC) simulations of DNA/RNA
helices using the base-pair level model, coded with Python. HelixMC is
distributed under the GPLv3 licence.

The project is authored by Fang-Chieh Chou in 2013, under the supervison of
Dr. Rhiju Das, at the Biochemistry Department of Stanford Unviersity.

Links
=====
- Source code: https://github.com/fcchou/HelixMC
- HTML documentation: http://helixmc.readthedocs.org/
- Das Lab @ Stanford: http://daslab.stanford.edu

Dependencies
============
The required dependencies to build the software are Python >= 2.7,
Numpy >= 1.6, Matplotlib >= 1.1.0, and a working C/C++ compiler.

Install
=======
The easiest way to install is to use pip install::

    $ pip install helixmc

Alternatively, one can download the source code from the latest GitHub
repository. And install using `setup.py`::

    $ python setup.py build
    $ sudo python setup.py install

Instead of installing using setup.py, you can just add your HelixMC folder
into the system's ``$PATH`` and ``$PYTHONPATH``. In bash this can be done by
adding the following lines to your ``~/.bashrc``::

    export PATH=$PATH:<HelixMC Path>
    export PYTHONPATH=$PYTHONPATH:<HelixMC Path>

Then build the Cython extension. Under the ``helixmc/`` folder, run::

    $ python _cython_build.py build_ext --inplace

Note that this requires you to have Cython installed. Otherwise you can choose
to build the c source file, then you do not need Cython::

    $ python _c_build.py build_ext --inplace
