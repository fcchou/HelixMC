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

import numpy as np
import abc
import os.path
from util import params2coords, _circmean
from __init__ import random


#####Random base-pair step parameters generator classes#####
class RandomStepBase(object):
    """
    Base class for Random bp-step generator, for inheritence only.

    See Also
    --------
    RandomStepSimple :
        Simple random bp-step generator.
    RandomStepAgg :
        Random bp-step generator by aggregrating multiple independent bp-step
        generators.
    """
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def __init__(self):
        return

    @abc.abstractmethod
    def __call__(self):
        return

    @abc.abstractproperty
    def params_avg(self):
        return

    @abc.abstractproperty
    def gaussian_sampling(self):
        return


class RandomStepSimple(RandomStepBase):
    '''
    Simple random bp-step generator.
    Pick random datasets from database or sample from a multivariate Gaussian
    built from the database.

    Parameters
    ----------
    params : ndarray, shape (N,6), optional
        Input base-pair step parameters. Usually obtained by analyzing
        the structures in PDB.
        Must be specified if params_cov and params_avg is not being input,
        or if gaussian_sampling is set to False.
        Order = [Shift, Slide, Rise, Tilt, Roll, Twist]
        Distance in unit of Å, angle in unit of radians.
    params_cov : ndarray, shape (6,6), optional
        Covariance matrix of the multivariate Gaussian for step parameters.
        Used for Gaussian sampling.
    params_avg : ndarray, shape (6,6), optional
        Averages of the multivariate Gaussian for step parameters.
        Used for Gaussian sampling.
    gaussian_sampling : bool, optional
        Whether to sample assuming a multivariate Gaussian. Default is True.

    Raises
    ------
    ValueError
        If params is not specified, and either params_avg and params_cov
        are not specified.
        If params is not specified, but gaussian_sampling is set to False.

    Attributes
    ----------
    `gaussian_sampling` : bool
        See the Parameters section.
    `params_avg` : ndarray
        See the Parameters section.
    `params_cov` : ndarray
        See the Parameters section.
    `params` : ndarray
        See the Parameters section. Return None if not specified at init.

    See Also
    --------
    RandomStepBase :
        Base class for Random bp-step generator, for inheritence only.
    RandomStepAgg :
        Random bp-step generator by aggregrating multiple independent bp-step
        generators.
    '''
    def __init__(
        self, params=None, params_cov=None,
        params_avg=None, gaussian_sampling=True
    ):
        self._gaussian_sampling = gaussian_sampling
        self.params = params
        if params is None:
            if params_avg is None or params_cov is None:
                raise ValueError(
                    'params is not specified; params_avg and'
                    ' params_cov are not specified.')
            self._params_cov = params_cov
            self._params_avg = params_avg
            ### To avoid float-point error ###
            self._params_avg[np.abs(self._params_avg) < 1e-15] = 0
            self._params_cov[np.abs(self._params_cov) < 1e-15] = 0

    @classmethod
    def load_gaussian_params(cls, filename):
        '''
        Load a single Gaussian parameter file from disk.

        Parameters
        ----------
        filename : str
            Name of the Gaussian parameter file, in .npy or plain text format.
            First row is the mean parameters and the following six rows
            represents the covariance matrix.
            The routine will look for helixmc/data folder if the file cannot
            be found.
            Order = [Shift, Slide, Rise, Tilt, Roll, Twist]
            Distance in unit of Å, angle in unit of radians.

        Returns
        -------
        random_step : RandomStepSimple
            Random step generator corresponds to the input file.

        Raises
        ------
        ValueError : if the input file does not exist.
        '''
        file_real = filename[:]
        if not os.path.exists(file_real):
            location = os.path.abspath(os.path.dirname(__file__))
            file_real = os.path.join(location, 'data/', file_real)
        if not os.path.exists(file_real):
            raise ValueError('Cannot find the input file %s!' % filename)
        if file_real[-3:] == 'npy':
            gaussian_params = np.load(file_real)
        else:
            gaussian_params = np.loadtxt(file_real)
        params_avg = gaussian_params[0]
        params_cov = gaussian_params[1:]
        return cls(
            params_avg=params_avg, params_cov=params_cov)

    @classmethod
    def load_params(cls, filename, gaussian_sampling=True):
        '''
        Load a simple bp-step parameter file from disk.

        Parameters
        ----------
        filename : str
            Name of the parameter file, in .npy or text format. It should
            contain a (N * 6) numpy array, where each row is a step parameter.
            The routine will look for helixmc/data folder if the file cannot
            be found.
            Order = [Shift, Slide, Rise, Tilt, Roll, Twist]
            Distance in unit of Å, angle in unit of radians.
        gaussian_sampling : bool, optional
            Whether to sample assuming a multivariate Gaussian.

        Returns
        -------
        random_step : RandomStepSimple
            Random step generator corresponds to the input file.

        Raises
        ------
        ValueError : if the input file does not exist.
        '''
        file_real = filename[:]
        if not os.path.exists(file_real):
            location = os.path.abspath(os.path.dirname(__file__))
            file_real = os.path.join(location, 'data/', file_real)
        if not os.path.exists(file_real):
            raise ValueError('Cannot find the input file %s!' % filename)
        if file_real[-3:] == 'npy':
            params = np.load(file_real)
        else:
            params = np.loadtxt(file_real)
        return cls(params, gaussian_sampling=gaussian_sampling)

    @property
    def gaussian_sampling(self):
        return self._gaussian_sampling

    @gaussian_sampling.setter
    def gaussian_sampling(self, val):
        self._gaussian_sampling = val
        if val:
            self._max_cached_data = 20000
            self._idx = self._max_cached_data
        elif self._params is None:
            raise ValueError(
                'Params is not specified, '
                'but gaussian_sampling is set to False.')

    @property
    def params_avg(self):
        return self._params_avg.copy()

    @property
    def params_cov(self):
        return self._params_cov.copy()

    @property
    def params(self):
        return self._params.copy()

    @params.setter
    def params(self, val):
        self._params = val
        if val is not None:
            self._params_avg = np.hstack((
                np.average(val[:, :3], axis=0),
                _circmean(val[:, 3:], axis=0)))
            val = val - self._params_avg
            val[:, 3:][val[:, 3:] > np.pi] -= 2 * np.pi
            val[:, 3:][val[:, 3:] <= -np.pi] += 2 * np.pi
            self._params_cov = np.cov(val, rowvar=0)
            ### To avoid float-point error ###
            self._params_avg[np.abs(self._params_avg) < 1e-15] = 0
            self._params_cov[np.abs(self._params_cov) < 1e-15] = 0
            ######
            self._o_list, self._R_list = params2coords(self._params)
            self._n_bp_step = self._params.shape[0]
        self.gaussian_sampling = self.gaussian_sampling

    def __call__(self):
        '''
        Draw one random bp-step parameter set from the distribution.

        Returns
        -------
        params: ndarray, shape (6)
            Randomly generated bp-step parameter set.
        o : ndarray, shape (3)
            Corresponding bp-center for the 2nd bp of the bp-step
            (1st bp is at origin and oriented with the coordinate frame).
        R : ndarray, shape (3,3)
            Corresponding frame for the 2nd bp of the bp-step.
        '''
        if self._gaussian_sampling:
            self._idx += 1
            if self._idx >= self._max_cached_data:
                self._idx = 0
                self._cached_params = random.multivariate_normal(
                    self._params_avg, self._params_cov, self._max_cached_data)
                self._o_list, self._R_list = params2coords(self._cached_params)
            return (
                self._cached_params[self._idx], self._o_list[self._idx],
                self._R_list[self._idx])
        else:
            i = random.randint(self._n_bp_step)
            return self._params[i], self._o_list[i], self._R_list[i]


class RandomStepAgg(RandomStepBase):
    '''
    Random bp-step generator by aggregating multiple
    independent simple bp-step generators.
    Useful for sequence dependence simulations.

    Parameters
    ----------
    data_file : str, optional
        Pre-curated database file with sequence dependence in .npz format.
    gaussian_sampling : bool, optional
        Whether to sample assuming a multivariate Gaussian. Default is True.

    Attributes
    ----------
    `gaussian_sampling` : bool
        See the Parameters section.
    `params_avg` : ndarray
        Average values for the step parameters distribution.
    `rand_list` : list
        List of all RandomStep objects in the aggregation.
    `names` : list
        List of all names of RandomStep in the aggregation.

    See Also
    --------
    RandomStepBase :
        Base class for Random bp-step generator, for inheritence only.
    RandomStepSimple :
        Simple random bp-step generator.
    '''
    def __init__(self, data_file=None, gaussian_sampling=True):
        self._names = []
        self._rand_list = []
        self._gaussian_sampling = gaussian_sampling
        if data_file is not None:
            self.load_from_file(data_file)

    def load_from_file(self, data_file):
        '''
        Load data file in .npz format and append to the current aggregation.
        The routine will look for helixmc/data folder if
        the data_file cannot be found.

        Parameters
        ----------
        data_file : str, optional
            Pre-curated database file with sequence dependence in .npz format.

        Raises
        ------
        ValueError : if the data_file does not exist.
        '''
        data_file_real = data_file[:]
        if not os.path.exists(data_file_real):
            location = os.path.abspath(os.path.dirname(__file__))
            data_file_real = os.path.join(location, 'data/', data_file_real)
        if not os.path.exists(data_file_real):
            raise ValueError('Cannot find the data file %s!' % data_file)
        data = np.load(data_file_real)
        for name in data.files:
            self.append_random_step(name, RandomStepSimple(
                params=data[name], gaussian_sampling=self.gaussian_sampling))

    def append_random_step(self, name, random_step):
        '''
        Append one additional random bp-step generator to the aggregation.

        Parameters
        ----------
        name : str
            Name of the random bp-step generator.
        random_step : subclass of RandomStepBase
            The random bp-step generator to be added to the aggregation

        Raises
        ------
        TypeError
            If random_step does not belong to a subclass of `RandomStepBase`.
        ValuError
            If the input name already exists in self.names.
        '''
        if not isinstance(random_step, RandomStepBase):
            raise TypeError(
                'Appended random_step is not a subclass of RandomStepBase.')
        if name in self._names:
            raise ValueError(
                "Appended random_step name '%s' already "
                "exists in list of names!" % name)
        self._names.append(name)
        self._rand_list.append(random_step)

    def symmetrize(self):
        '''
        Symmetrize the dataset by counting from opposite direction.

        Raises
        ------
        ValueError :
            If the dataset contain mixed RNA and DNA bases (U and T).
        '''
        is_DNA = False
        is_RNA = False
        for name in self._names:
            if 'T' in name:
                is_DNA = True
            if 'U' in name:
                is_RNA = True
        if is_DNA and is_RNA:
            raise ValueError(
                "The data is RNA-DNA mixture!! Cannot symmetrize.")
        if is_DNA:
            name_conv = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        else:
            name_conv = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}

        new_params = {}
        for name in self._names:
            sym_name = ''
            for lett in name:
                sym_name += name_conv[lett]
            sym_name = sym_name[::-1]
            params1 = self.get_rand_step(name).params
            params1[:, 0] *= -1
            params1[:, 3] *= -1
            params2 = np.vstack((self.get_rand_step(sym_name).params, params1))
            new_params[sym_name] = params2

        for name in self._names:
            self.get_rand_step(name).params = new_params[name]

    def clear_all(self):
        'Clear all random bp-step generator in the aggregation.'
        self._names[:] = []
        self._rand_list[:] = []

    @property
    def gaussian_sampling(self):
        return self._gaussian_sampling

    @gaussian_sampling.setter
    def gaussian_sampling(self, val):
        self._gaussian_sampling = val
        for i in self._rand_list:
            i.gaussian_sampling = val

    @property
    def params_avg(self):
        params_list = np.empty((len(self._rand_list), 6))
        for i, rand in enumerate(self._rand_list):
            params_list[i] = rand.params_avg
        return np.hstack((
            np.average(params_list[:, :3], axis=0),
            _circmean(params_list[:, 3:], axis=0)))

    @property
    def names(self):
        return self._names[:]

    @property
    def rand_list(self):
        return self.rand_list[:]

    def name2rand_idx(self, name):
        '''
        Get the index of a RandomStep object in rand_list from its name.

        Parameters
        ----------
        name : str
            Name of the RandomStep in the aggregation.

        Returns
        -------
        idx : int
            The corresponding index of the RandomStep object.
        '''
        return self._names.index(name)

    def get_rand_step(self, identifier):
        '''
        Return one RandomStep object stored in the aggregation.

        Parameters
        ----------
        identifier : int or str
            Index or name of the requested random step generator.

        Returns
        -------
        rand_step : subclass of RandomStepBase
             RandomStep object stored in the aggregation.

        Raises
        ------
        TypeError :
            If the input identifier is not int or str
        '''
        if isinstance(identifier, int):
            return self._rand_list[identifier]
        elif isinstance(identifier, str):
            return self._rand_list[self.name2rand_idx(identifier)]
        else:
            raise TypeError(
                "Input identifier has invalid datatype "
                "(is %s, should be int or str)!!" % type(identifier))

    def __call__(self, identifier=None):
        '''
        Draw one random bp-step parameter set from the distribution.
        Randomly select one generator in the aggregation and
        return its generated result if no input parameter is given.

        Parameters
        ----------
        identifier : int or str
            Index or name of the requested random step generator.

        Returns
        -------
        params: ndarray, shape (6)
            Randomly generated bp-step parameter set.
        o : ndarray, shape (3)
            Corresponding bp-center for the 2nd bp of the bp-step (1st bp is
            at origin and oriented with the coordinate frame).
        R : ndarray, shape (3,3)
            Corresponding frame for the 2nd bp of the bp-step.

        Raises
        ------
        TypeError :
            If the input identifier is not int or str
        '''
        if identifier is None:
            return self._rand_list[random.randint(len(self._rand_list))]()
        elif isinstance(identifier, int):
            return self._rand_list[identifier]()
        elif isinstance(identifier, str):
            return self._rand_list[self.name2rand_idx(identifier)]()
        else:
            raise TypeError(
                "Input identifier has invalid datatype "
                "(is %s, should be int or str)!!" % type(identifier))
