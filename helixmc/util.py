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
import math
from numpy import pi
from . import random, kBT
from os.path import exists
from _util_cython import writhe_fuller, writhe_exact, ribbon_twist


#####Utility functions#####
def Rz(theta):
    '''
    Return z-rotation matrices with rotational angle theta.

    Parameters
    ----------
    theta : array-like
        Rotation angles of the matrix in radians.

    Returns
    -------
    rot_matrix : ndarray
        Corresponding z-rotation martrices for each input angle,
        align with the first index.

    See Also
    --------
    Ry : Return y-rotation matrices with rotational angle theta.
    Rx : Return x-rotation matrices with rotational angle theta.
    R_axis : Return rotation matrices with rotational angle theta
    along an arbitary rotation axis.
    '''
    if isinstance(theta, np.ndarray):  # array version
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        rot_matrix = np.zeros((theta.shape[0], 3, 3))
        rot_matrix[:, 0, 0] = cos_theta
        rot_matrix[:, 0, 1] = -sin_theta
        rot_matrix[:, 1, 0] = sin_theta
        rot_matrix[:, 1, 1] = cos_theta
        rot_matrix[:, 2, 2] = 1
    else:  # Scalar version
        sin_theta = math.sin(theta)
        cos_theta = math.cos(theta)
        rot_matrix = np.eye(3)
        rot_matrix[0, 0] = cos_theta
        rot_matrix[0, 1] = -sin_theta
        rot_matrix[1, 0] = sin_theta
        rot_matrix[1, 1] = cos_theta
    return rot_matrix


def Rx(theta):
    '''
    Return x-rotation matrices with rotational angle theta.

    Please refer to the documentation for `Rz` for further details.

    See Also
    --------
    Rz, Ry, R_axis
    '''
    if isinstance(theta, np.ndarray):  # array version
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        rot_matrix = np.zeros((theta.shape[0], 3, 3))
        rot_matrix[:, 1, 1] = cos_theta
        rot_matrix[:, 1, 2] = -sin_theta
        rot_matrix[:, 2, 1] = sin_theta
        rot_matrix[:, 2, 2] = cos_theta
        rot_matrix[:, 0, 0] = 1
    else:  # Scalar version
        sin_theta = math.sin(theta)
        cos_theta = math.cos(theta)
        rot_matrix = np.eye(3)
        rot_matrix[1, 1] = cos_theta
        rot_matrix[1, 2] = -sin_theta
        rot_matrix[2, 1] = sin_theta
        rot_matrix[2, 2] = cos_theta
    return rot_matrix


def Ry(theta):
    '''
    Return y-rotation matrices with rotational angle theta.

    Please refer to the documentation for `Rz` for further details.

    See Also
    --------
    Rz, Rx, R_axis
    '''
    if isinstance(theta, np.ndarray):  # array version
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        rot_matrix = np.zeros((theta.shape[0], 3, 3))
        rot_matrix[:, 0, 0] = cos_theta
        rot_matrix[:, 0, 2] = sin_theta
        rot_matrix[:, 2, 0] = -sin_theta
        rot_matrix[:, 2, 2] = cos_theta
        rot_matrix[:, 1, 1] = 1
    else:  # Scalar version
        sin_theta = math.sin(theta)
        cos_theta = math.cos(theta)
        rot_matrix = np.eye(3)
        rot_matrix[0, 0] = cos_theta
        rot_matrix[0, 2] = sin_theta
        rot_matrix[2, 0] = -sin_theta
        rot_matrix[2, 2] = cos_theta
    return rot_matrix


def R_axis(theta, axis):
    '''
    Return rotation matrices with rotational angle theta along
    an arbitary rotation axis.

    Parameters
    ----------
    theta : array-like
        Rotation angles of the matrix in radians.
    axis : ndarray, shape (N,3)
        Rotational axis for the rotation being performed, align with first
        index (i.e. axis[3] is the rotational axis for angle theta[3]).

    Returns
    -------
    rot_matrix : ndarray
        Corresponding rotation martrices for each input angle,
        align with the first index.

    See Also
    --------
    Rz :
        Return z-rotation matrices with rotational angle theta.
    Ry :
        Return y-rotation matrices with rotational angle theta.
    Rx :
        Return x-rotation matrices with rotational angle theta.
    '''
    if isinstance(theta, np.ndarray):  # array version
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)
        axis = axis / np.sqrt(np.sum(axis * axis, axis=1))[:, np.newaxis]
        rot_matrix = np.zeros((theta.shape[0], 3, 3))
        rot_matrix[:, 0, 0] = cos_theta + axis[:, 0] ** 2 * (1 - cos_theta)
        rot_matrix[:, 1, 1] = cos_theta + axis[:, 1] ** 2 * (1 - cos_theta)
        rot_matrix[:, 2, 2] = cos_theta + axis[:, 2] ** 2 * (1 - cos_theta)
        rot_matrix[:, 0, 1] = (
            axis[:, 0] * axis[:, 1] * (1 - cos_theta) -
            axis[:, 2] * sin_theta)
        rot_matrix[:, 1, 0] = (
            axis[:, 0] * axis[:, 1] * (1 - cos_theta) +
            axis[:, 2] * sin_theta)
        rot_matrix[:, 0, 2] = (
            axis[:, 0] * axis[:, 2] * (1 - cos_theta) +
            axis[:, 1] * sin_theta)
        rot_matrix[:, 2, 0] = (
            axis[:, 0] * axis[:, 2] * (1 - cos_theta) -
            axis[:, 1] * sin_theta)
        rot_matrix[:, 1, 2] = (
            axis[:, 1] * axis[:, 2] * (1 - cos_theta) -
            axis[:, 0] * sin_theta)
        rot_matrix[:, 2, 1] = (
            axis[:, 1] * axis[:, 2] * (1 - cos_theta) +
            axis[:, 0] * sin_theta)
    else:  # scalar
        sin_theta = math.sin(theta)
        cos_theta = math.cos(theta)
        axis = axis / math.sqrt(np.dot(axis, axis))
        rot_matrix = np.identity(3)
        rot_matrix[0, 0] = cos_theta + axis[0] ** 2 * (1 - cos_theta)
        rot_matrix[1, 1] = cos_theta + axis[1] ** 2 * (1 - cos_theta)
        rot_matrix[2, 2] = cos_theta + axis[2] ** 2 * (1 - cos_theta)
        rot_matrix[0, 1] = (
            axis[0] * axis[1] * (1 - cos_theta) -
            axis[2] * sin_theta)
        rot_matrix[1, 0] = (
            axis[0] * axis[1] * (1 - cos_theta) +
            axis[2] * sin_theta)
        rot_matrix[0, 2] = (
            axis[0] * axis[2] * (1 - cos_theta) +
            axis[1] * sin_theta)
        rot_matrix[2, 0] = (
            axis[0] * axis[2] * (1 - cos_theta) -
            axis[1] * sin_theta)
        rot_matrix[1, 2] = (
            axis[1] * axis[2] * (1 - cos_theta) -
            axis[0] * sin_theta)
        rot_matrix[2, 1] = (
            axis[1] * axis[2] * (1 - cos_theta) +
            axis[0] * sin_theta)
    return rot_matrix


def params2coords(params):
    '''
    Convert base-pair step parameters to bp-center coordinates
    and rotation matrix.

    Parameters
    ----------
    params : ndarray, shape (N,6)
        Input base-pair step parameters. Distances in unit of Å,
        angles in unit of radians.
        Order = [Shift, Slide, Rise, Tilt, Roll, Twist]

    Returns
    -------
    o2 : ndarray, shape (N,3)
        Coordiantes for the bp-centers of the 2nd base-pair.
    R2 : ndarray, shape (N,3,3)
        Rotational matrix for the transformation.

    See Also
    --------
    coords2params :
        Convert bp-center coordinates and rotation matrix to base-pair step
        parameters.
    '''
    if len(params.shape) == 1:  # 1D case
        Tau = math.sqrt(params[3] ** 2 + params[4] ** 2)
        Phi = math.atan2(params[3], params[4])
        o2 = (
            Rz(params[5] * 0.5 - Phi).dot(Ry(Tau * 0.5)).
            dot(Rz(Phi)).dot(params[0:3]))
        R2 = (
            Rz(-Phi + params[5] * 0.5).dot(Ry(Tau)).
            dot(Rz(params[5] * 0.5 + Phi)))
    else:  # 2D case
        Tau = np.sqrt(params[:, 3] ** 2 + params[:, 4] ** 2)
        Phi = np.arctan2(params[:, 3], params[:, 4])
        o2 = np.einsum(
            'ijk,ikl,ilm,im->ij', Rz(params[:, 5] * 0.5 - Phi),
            Ry(Tau * 0.5), Rz(Phi), params[:, 0:3])
        R2 = np.einsum(
            'ijk,ikl,ilm->ijm', Rz(-Phi + params[:, 5] * 0.5),
            Ry(Tau), Rz(params[:, 5] * 0.5 + Phi))
    return o2, R2


def coords2params(o2, R2):
    '''
    Convert bp-center coordinates and rotation matrix to
    base-pair step parameters.
    The frame of bp1 is used as reference ( o1 = [0 0 0] and R1 = np.eye(3) )

    Parameters
    ----------
    o2 : ndarray, shape (N,3)
        Coordiantes for the bp-centers of the 2nd base-pair.
    R2 : ndarray, shape (N,3,3)
        Rotational matrix for the transformation.

    Returns
    -------
        Base-pair step parameters.
        Distance in unit of Å, angle in unit of radians.
        Order = [Shift, Slide, Rise, Tilt, Roll, Twist]

    See Also
    --------
    params2coords :
        Convert base-pair step parameters to bp-center coordinates and
        rotation matrix.
    '''
    if len(o2.shape) == 1:  # 1D case
        Tau = math.acos(R2[2, 2])
        if math.sin(Tau) > 0:
            alpha = math.atan2(R2[1, 2], R2[0, 2])
            gamma = math.atan2(R2[2, 1], -R2[2, 0])
        else:
            alpha = math.atan2(-R2[1, 2], -R2[0, 2])
            gamma = math.atan2(-R2[2, 1], R2[2, 0])
        Twst = alpha + gamma
        Phi = 0.5 * (gamma - alpha)
        if Twst > pi:
            Twst -= 2 * pi
            Phi -= pi
        elif Twst <= -pi:
            Twst += 2 * pi
            Phi += pi
        Roll = Tau * np.cos(Phi)
        Tilt = Tau * np.sin(Phi)
        if Tau < 1e-6:
            Twst = math.atan2(R2[1, 0], R2[0, 0])
        d_vec = Rz(-Phi).dot(Ry(-Tau * 0.5)).dot(Rz(-Twst * 0.5 + Phi)).dot(o2)
        return np.hstack((d_vec, Tilt, Roll, Twst))
    else:  # 2D case
        Tau = np.arccos(R2[:, 2, 2])
        sign_sin_Tau = np.sign(np.sin(Tau))
        alpha = np.arctan2(
            sign_sin_Tau * R2[:, 1, 2], sign_sin_Tau * R2[:, 0, 2])
        gamma = np.arctan2(
            sign_sin_Tau * R2[:, 2, 1], -sign_sin_Tau * R2[:, 2, 0])
        Twst = alpha + gamma
        Phi = 0.5 * (gamma - alpha)
        Phi[Twst > pi] += pi
        Phi[Twst <= -pi] += pi
        Twst[Twst > pi] -= 2 * pi
        Twst[Twst <= -pi] += 2 * pi
        Roll = Tau * np.cos(Phi)
        Tilt = Tau * np.sin(Phi)
        special_idx = (Tau < 1e-6)
        if np.any(special_idx):
            Twst[special_idx] = np.arctan2(
                R2[special_idx, 1, 0], R2[special_idx, 0, 0])
        d_vec = np.einsum(
            'ijk,ikl,ilm,im->ij', Rz(-Phi), Ry(-Tau * 0.5),
            Rz(-Twst * 0.5 + Phi), o2)
        return np.column_stack((d_vec, Tilt, Roll, Twst))


def coords2dr(coord):
    '''
    Convert xyz coordinates of axis curve to delta-r vectors.

    Parameters
    ----------
    coord : ndarray, shape (N,3)
        The axis curve xyz coordinates in unit of Å (r vectors).

    Returns
    -------
    dr : ndarray, shape (N,3)
        Delta-r vectors (dr[i] = r[i+1] - r[i]).

    See Also
    --------
    dr2coords :
        Convert delta-r vectors to coordinates.
    '''
    dr = coord[1:] - coord[:-1]
    return dr


def dr2coords(dr):
    '''
    Convert delta-r vectors to coordinates.

    Parameters
    ----------
    dr : ndarray, shape (N,3)
        Input delta_r vectors (dr[i] = r[i+1] - r[i]).

    Returns
    -------
    coord : ndarray, shape (N+1,3)
        Coordinates for the base-pair-centers (r[i]).

    See Also
    --------
    coords2dr :
        Convert xyz coordinates of axis curve to delta-r vectors.
    '''
    coord = np.vstack((np.zeros(3), np.cumsum(dr, axis=0)))
    return coord


def params2data(params, frame0=None):
    '''
    Convert step parameters to delta-r vectors and frames.

    Parameters
    ----------
    params : ndarray, shape (N,6)
        Input base-pair step parameters.
        Order = [Shift, Slide, Rise, Tilt, Roll, Twist]
        Distance in unit of Å, angle in unit of radians.
    frame0 : ndarray, shape (3,3), optional
        The frame for the 1st base-pair. Default set to np.eye(3).

    Returns
    -------
    dr : ndarray, shape (N,3)
        Delta-r vectors (dr[i] = r[i+1] - r[i]).
    frames : ndarray, shape (N+1,3,3)
        Frame of each base-pair.

    See Also
    --------
    data2params :
        Convert delta-r vectors and frames to step parameters.
    '''
    if frame0 is None:
        frame0 = np.eye(3)
    n_bp = params.shape[0] + 1
    frames = np.zeros((n_bp, 3, 3))
    o, R = params2coords(params)
    unitarize(R)
    for i in xrange(n_bp):
        if i == 0:
            frames[i] = unitarize(frame0)
        else:
            frames[i] = unitarize(frames[i-1].dot(R[i-1]))
    dr = np.einsum('ijk,ik->ij', frames[:-1], o)
    return dr, frames


def _data2params(dr, f1, f2):
    '''
    Convert delta-r, frame for 1st and 2nd bp to step parameters.
    '''
    if len(dr.shape) == 1:  # 1D case
        return coords2params(f1.T.dot(dr), f1.T.dot(f2))
    else:
        o = np.einsum('ijk,ij ->ik', f1, dr)
        R = np.einsum('ikj,ikl->ijl', f1, f2)
        params = coords2params(o, unitarize(R))
        return params


def data2params(dr, frames):
    '''
    Convert delta-r vectors and frames to step parameters.

    Parameters
    ----------
    dr : ndarray, shape (N,3)
        Delta-r vectors (dr[i] = r[i+1] - r[i]).
    frames : ndarray, shape (N+1,3,3)
        Frame of each base-pair.

    Returns
    -------
    params : ndarray, shape (N,6)
        Base-pair step parameters.
        Distance in unit of Å, angle in unit of radians.
        Order = [Shift, Slide, Rise, Tilt, Roll, Twist]

    See Also
    --------
    params2data :
        Convert step parameters to delta-r vectors and frames.
    '''
    return _data2params(dr, frames[:-1], frames[1:])


def frames2params(o1, o2, f1, f2):
    '''
    Convert bp coordinates and frames to step parameters.

    Parameters
    ----------
    o1 : ndarray, shape (N,3)
        Origins for the bp-centers of the 1st base-pair.
    o2 : ndarray, shape (N,3)
        Origins for the bp-centers of the 2nd base-pair.
    f1 : ndarray, shape (N,3,3)
        Frames for the 1st base-pair.
    f2 : ndarray, shape (N,3,3)
        Frames for the 2nd base-pair.

    Returns
    -------
    params : ndarray, shape (N,6)
        Base-pair step parameters.
        Distance in unit of Å, angle in unit of radians.
        Order = [Shift, Slide, Rise, Tilt, Roll, Twist]

    See Also
    --------
    frames2params_3dna :
        Convert bp coordinates and frames in 3DNA format to step parameters.
    '''
    return _data2params(o2 - o1, f1, f2)


def frames2params_3dna(o1, o2, f1, f2):
    '''
    Convert bp coordinates and frames in 3DNA format to step parameters.
    Note that the 3DNA base-pair frames differ from HelixMC
    by a transpose operation.

    See Also
    --------
    frames2params :
        Convert bp coordinates and frames to step parameters.
    '''
    if len(o1.shape == 1):  # 1D case
        return frames2params(o1, o2, f1.T, f2.T)
    else:
        return frames2params(
            o1, o2, np.transpose(f1, (0, 2, 1)),
            np.transpose(f2, (0, 2, 1)))


def params_join(params):
    '''
    Convert consequtive bp-params to the params between 1st and last bp.

    Parameters
    ----------
    params : ndarray, shape (N,6)
        Input base-pair step parameters. Distances in unit of Å,
        angles in unit of radians.
        Order = [Shift, Slide, Rise, Tilt, Roll, Twist]

    Returns
    -------
    params :  ndarray, shape (6)
        Base-pair step parameters between 1st and last.
    '''
    if len(params.shape) == 1:
        return params
    if params.shape[0] == 1:
        return params[0]
    dr, frames = params2data(params)
    o = np.sum(dr, axis=0)
    R = frames[-1]
    return coords2params(o, R)


def unitarize(R):
    '''
    Enforce unitarity of the input matrix using Gram-Schmidt process.
    This function modifies the input matrix inplace.

    Parameters
    ----------
    R : ndarray, shape (N,3,3)
        (List of) input matrix to be unitarized.

    Returns
    -------
    R : ndarray, shape (N,3,3)
        Unitarized input matrix.
    '''
    if len(R.shape) == 2:  # 1D case
        R[0] /= math.sqrt(R[0].dot(R[0]))
        R[1] -= R[1].dot(R[0]) * R[0]
        R[1] /= math.sqrt(R[1].dot(R[1]))
        R[2] -= R[2].dot(R[0]) * R[0]
        R[2] -= R[2].dot(R[1]) * R[1]
        R[2] /= math.sqrt(R[2].dot(R[2]))
    else:  # 2D case
        R[:, 0] /= np.sqrt(np.sum(R[:, 0] ** 2, axis=1))[:, np.newaxis]
        R[:, 1] -= (
            np.einsum('ij,ij->i', R[:, 1], R[:, 0])[:, np.newaxis] * R[:, 0])
        R[:, 1] /= np.sqrt(np.sum(R[:, 1] ** 2, axis=1))[:, np.newaxis]
        R[:, 2] -= (
            np.einsum('ij,ij->i', R[:, 2], R[:, 0])[:, np.newaxis] * R[:, 0])
        R[:, 2] -= (
            np.einsum('ij,ij->i', R[:, 2], R[:, 1])[:, np.newaxis] * R[:, 1])
        R[:, 2] /= np.sqrt(np.sum(R[:, 2] ** 2, axis=1))[:, np.newaxis]
    return R


def MC_acpt_rej(score_old, score_new, kT=kBT):
    '''
    Decide whether to accept a Monte-Carlo move.

    Parameters
    ----------
    score_old : float
        score before the proposed update.
    score_new : float
        score after the proposed update.
    kT : float
        Temperature times Boltzmann constant, in pN.Å.

    Returns
    -------
    out : bool
        True for accept, False for reject.
    '''
    exponent = (score_new - score_old) / kT
    if exponent <= 1e-9 or math.exp(-exponent) >= random.random_sample():
        return True
    else:
        return False


def read_seq_from_fasta(fasta):
    '''
    Read the sequence from a fasta file.

    Parameters
    ----------
    fasta : str
        Name of the fasta file.

    Returns
    -------
    seq : str
        Sequence stored in fasta.

    Raises
    ------
    ValueError
        If the fasta file does not exist.
    '''
    if not exists(fasta):
        raise ValueError('The fasta file %s does not exist.' % fasta)
    comment_symbols = ['>', ';', '#']
    seq = ''
    for line in open(fasta):
        if line[0] not in comment_symbols:
            seq += line.strip()
    return seq


def _circmean(arr, axis=None):
    '''
    Circular mean of angles.
    Results are in [-pi, pi]. Note that scipy.stats has a similar function,
    but we re-implement it because the scipy one depends on its version.

    Parameters
    ----------
    arr : ndarray
        Input numpy array.

    axis : int, optional
        Axis on which the mean is computed.

    Returns
    -------
    mean : ndarray
        The circular mean.
    '''
    sin_mean = np.average(np.sin(arr), axis=axis)
    cos_mean = np.average(np.cos(arr), axis=axis)
    return np.arctan2(sin_mean, cos_mean)
