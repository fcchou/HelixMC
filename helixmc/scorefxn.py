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

import abc
import warnings

#####Score function#####
class ScorefxnBase(object):
    '''
    Base class for scoring fucntion, for inheritence only.

    See Also
    --------
    ScorefxnTweezers : Score function for tweezers experiments.
    '''
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __init__(self):
        return

    @abc.abstractmethod
    def __call__(self, pose):
        return

class ScorefxnTweezers(ScorefxnBase):
    '''
    Score function for tweezers experiments.

    Parameters
    ----------
    force : float, optional
        Applied z-direction force to the helix, in pN.
    trap_stiffness : float, optional
        The stiffness of the angular trap, in pN.A.
    target_link : float, optional
        Center of the harmonic angular trap (linking number trap), in radians.

    Attributes
    ----------
    `force` : float
    `trap_stiffness` : float
    `target_link` : float
        See Parameters section above.

    See Also
    --------
    ScorefxnBase : Base class for scoring fucntion, for inheritence only.
    '''
    def __init__(self, force=0.0, trap_stiffness=0.0, target_link=None):
        self.force = force
        self.trap_stiffness = trap_stiffness
        self.target_link = target_link

    def __call__(self, pose):
        '''
        Score the input pose.

        Parameters
        ----------
        pose : HelixPose
            Input pose for scoring.

        Returns
        -------
        score : float
            Score of the pose.
        '''
        score = 0
        if self.force != 0:
            score -= pose.z_terminal * self.force
        if self.trap_stiffness != 0 and self.target_link != None:
            score += 0.5 * self.trap_stiffness * (pose.link_fuller - self.target_link) ** 2
            if not pose.compute_tw_wr:
                warnings.warn('pose.compute_tw_wr should be set to Ture for repeating scoring with target link!!!', RuntimeWarning)
        return score


