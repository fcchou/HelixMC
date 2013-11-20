#This file is part of HelixMC.
#    Copyright (C) 2013  Fang-Chieh Chou <fcchou@stanford.edu>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""
#######
HelixMC
#######

HelixMC is a software package for Monte-Carlo (MC) simulations of DNA/RNA
helices using the base-pair level model, coded with Python. HelixMC is
distributed under the GPLv3 licence.

The project is authored by Fang-Chieh Chou in 2013, under the supervison of
Dr. Rhiju Das, at the Biochemistry Department of Stanford Unviersity.

Contacts
--------
- Fang-Chieh Chou <fcchou@stanford.edu>
- Rhiju Das <rhiju@stanford.edu>

Links
-----
- Source code: https://github.com/fcchou/HelixMC
- HTML documentation: http://helixmc.readthedocs.org/
- Das Lab @ Stanford: http://daslab.stanford.edu
"""
#version
__version__ = '0.5'

#list of all the modules (files)
__all__ = [
    "util",
    "pose",
    "random_step",
    "score_function",
    "fit_function"
]

#Initialize Random number generator used throughout
import numpy as np
from numpy import random


def constant_seed(seed=24601):
    '''
    Set constant random seed.

    Parameters
    ----------
    seed : int, optional
        Random seed.
    '''
    random.seed(seed)
#####Useful constants#####
kB = 1.3806488e-1  # Boltzmann constant in pN.A/K
kBT = kB * 298.15  # kB.T at room temperature (25 degree Celsius)
ez = np.array([0, 0, 1])  # z unit vector
