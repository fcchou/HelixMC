# -*- coding: UTF-8 -*-

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
import numpy as np


#####Score function#####
class ScoreBase(object):
    '''
    Base class for scoring fucntion, for inheritence only.
    '''
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __init__(self):
        return

    @abc.abstractmethod
    def __call__(self, pose):
        return


class ScoreExt(ScoreBase):
    '''
    Score function for force-extension along Z-axis.

    Parameters
    ----------
    force : float
        Applied z-direction force to the helix, in pN.

    Attributes
    ----------
    `force` : float
        See Parameters section above.
    '''
    def __init__(self, force):
        self.force = force

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
        return -pose.z_terminal * self.force


class ScoreTorsionTrap(ScoreBase):
    '''
    Score function for torsional trap.

    Parameters
    ----------
    stiffness : float
        The stiffness of the torsional trap, in pN.Å.
    target_link : float
        Center of the harmonic torsional trap (link), in radians.

    Attributes
    ----------
    `stiffness` : float
    `target_link` : float
        See Parameters section above.
    '''
    def __init__(self, stiffness, target_link):
        self.stiffness = stiffness
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
        if not pose.compute_tw_wr:
            warnings.warn(
                'pose.compute_tw_wr should be set to Ture for repeating'
                'scoring with target link!!!', RuntimeWarning)
        return (
            0.5 * self.stiffness *
            (pose.link_fuller - self.target_link) ** 2)


class ScoreXyTrap(ScoreBase):
    '''
    Score function for xy trap.

    Parameters
    ----------
    stiffness : float
        The stiffness of the xy trap, in pN/Å.

    Attributes
    ----------
    `stiffness` : float
        See Parameters section above.
    '''
    def __init__(self, stiffness):
        self.stiffness = stiffness

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
        xy = pose.coord_terminal[:2]
        dist_sq = np.sum(xy ** 2)
        return 0.5 * self.stiffness * dist_sq


class ScoreAgg(ScoreBase):
    '''
    Score function aggregates of multiple score terms.

    Parameters
    ----------
    score_list : list, optional
        List of score terms (subclass of ScoreBase) in this score.

    Attributes
    ----------
    `score_list` : list
        See Parameters section above.
    `is_empty` : bool
        If the score_list is empty.
    '''
    def __init__(self, score_list=[]):
        self.score_list = score_list

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
        for term in self.score_list:
            score += term(pose)
        return score

    def append(self, score_term):
        '''
        Append new score term.

        Parameters
        ----------
        score_term : subclass of ScoreBase
            Score term to be appended.
        '''
        self.score_list.append(score_term)

    def clear(self):
        '''
        Clear the score_list.
        '''
        self.score_list = []

    @property
    def is_empty(self):
        return (not self.score_list)


class ScoreTweezers(ScoreAgg):
    '''
    Score function for tweezers experiments.

    Parameters
    ----------
    force : float, optional
        Applied z-stretching force, in pN.
    torsional_stiffness : float, optional
        Stiffness of the torsional trap, in pN.Å.
    target_link : float, optional
        Center of the torsional trap (link), in radians.
    xy_stiffness : float, optional
        Stiffness of xy trap, in pN/Å.
    '''
    def __init__(
            self,
            force=0,
            torsional_stiffness=0,
            target_link=None,
            xy_stiffness=0
    ):
        self.score_list = []
        if force != 0:
            self.score_list.append(ScoreExt(force))
        if torsional_stiffness != 0 and target_link is not None:
            self.score_list.append(
                ScoreTorsionTrap(torsional_stiffness, target_link))
        if xy_stiffness != 0:
            self.score_list.append(ScoreXyTrap(xy_stiffness))
