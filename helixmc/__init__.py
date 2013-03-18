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
HelixMC:

Author:
Fang-Chieh Chou <fcchou@stanford.edu>
"""

# list of all the modules (files)
__all__ = [
"util",
"pose",
"random_step",
"score_function"
]

#Initialize Random number generator used throughout
import numpy as np
from numpy import random

def constant_seed( seed=24601 ):
    '''
    Set constant random seed.

    Parameters
    ---------
    seed : int, optional
        Random seed.
    '''
    random.seed(seed)

#####Useful constants#####
kB = 1.3806488e-1 #Boltzmann constant in pN.A/K
kBT = kB * 298.15 #kB.T at room temperature (25 degree Celsius)
ez = np.array( [0,0,1] ) #z unit vector
from math import pi
