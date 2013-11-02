.. HelixMC documentation master file, created by
   sphinx-quickstart on Thu Mar 21 15:51:36 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

HelixMC documentation
=====================

HelixMC is a software package for Monte-Carlo (MC) simulations of DNA/RNA
helices using the base-pair level model. It provides a powerful tool to
understand the flexibility of DNA/RNA helices through numerical simualtions.

The base-pair level model, first developed by Olson and collegues, bridges
between the simple elastic rod model and the full-atom representation,
providing a reasonably sophiscated and yet computationally tractable way to
model long DNA/RNA helices up to thousands of base-pairs and to evalute
their mechanical properties. HelixMC has the utility of applying external
stetching forces and torques to the helix, and measuring the helix extension
and the rotation of the helix (known as the linking number) during the
process. This properly emulate the setup of recent single-molecule tweezers
experiments, making HelixMC a useful tool for direct simulations of these
experiments.

HelixMC is coded in Python in an object-oriented fashion. Rate-limiting core
computations are speeded up using Cython. Therefore HelixMC provides a
framework that is easy to use and to extend, as well as being reasonably fast
when it comes to large-scale computations.

.. toctree::
   :maxdepth: 2

   tutorial
   reference
   about
