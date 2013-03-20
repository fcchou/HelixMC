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
from numpy import pi
import math
import scipy.weave as weave
from __init__ import random, pi, kBT, ez
from os.path import exists

#####Utility functions#####
def Rz( theta ):
    '''
    Return z-rotation matrices with rotational angle theta.

    Parameters
    ----------
    theta : array-like
        Rotation angles of the matrix in radians.

    Returns
    -------
    rot_matrix : ndarray
        Corresponding z-rotation martrices for each input angle, align with the first index.

    See Also
    --------
    Ry : Return y-rotation matrices with rotational angle theta.
    Rx : Return x-rotation matrices with rotational angle theta.
    R_axis : Return rotation matrices with rotational angle theta along an arbitary rotation axis.
    '''
    if isinstance(theta, np.ndarray): # array version
        sin_theta = np.sin( theta )
        cos_theta = np.cos( theta )
        rot_matrix = np.zeros( (theta.shape[0], 3, 3) )
        rot_matrix[:,0,0] =  cos_theta
        rot_matrix[:,0,1] = -sin_theta
        rot_matrix[:,1,0] =  sin_theta
        rot_matrix[:,1,1] =  cos_theta
        rot_matrix[:,2,2] =  1
    else: #Scalar version
        sin_theta = math.sin( theta )
        cos_theta = math.cos( theta )
        rot_matrix = np.eye(3)
        rot_matrix[0,0] =  cos_theta
        rot_matrix[0,1] = -sin_theta
        rot_matrix[1,0] =  sin_theta
        rot_matrix[1,1] =  cos_theta
    return rot_matrix

def Rx( theta ):
    '''
    Return x-rotation matrices with rotational angle theta.

    Please refer to the documentation for `Rz` for further details.

    See Also
    --------
    Rz, Ry, R_axis
    '''
    if isinstance(theta, np.ndarray): # array version
        sin_theta = np.sin( theta )
        cos_theta = np.cos( theta )
        rot_matrix = np.zeros( (theta.shape[0], 3, 3) )
        rot_matrix[:,1,1] =  cos_theta
        rot_matrix[:,1,2] = -sin_theta
        rot_matrix[:,2,1] =  sin_theta
        rot_matrix[:,2,2] =  cos_theta
        rot_matrix[:,0,0] =  1
    else: #Scalar version
        sin_theta = math.sin( theta )
        cos_theta = math.cos( theta )
        rot_matrix = np.eye(3)
        rot_matrix[1,1] =  cos_theta
        rot_matrix[1,2] = -sin_theta
        rot_matrix[2,1] =  sin_theta
        rot_matrix[2,2] =  cos_theta
    return rot_matrix

def Ry( theta ):
    '''
    Return y-rotation matrices with rotational angle theta.

    Please refer to the documentation for `Rz` for further details.

    See Also
    --------
    Rz, Rx, R_axis
    '''
    if isinstance(theta, np.ndarray): # array version
        sin_theta = np.sin( theta )
        cos_theta = np.cos( theta )
        rot_matrix = np.zeros( (theta.shape[0], 3, 3) )
        rot_matrix[:,0,0] =  cos_theta
        rot_matrix[:,0,2] =  sin_theta
        rot_matrix[:,2,0] = -sin_theta
        rot_matrix[:,2,2] =  cos_theta
        rot_matrix[:,1,1] =  1
    else: #Scalar version
        sin_theta = math.sin( theta )
        cos_theta = math.cos( theta )
        rot_matrix = np.eye(3)
        rot_matrix[0,0] =  cos_theta
        rot_matrix[0,2] =  sin_theta
        rot_matrix[2,0] = -sin_theta
        rot_matrix[2,2] =  cos_theta
    return rot_matrix

def R_axis( theta, axis ):
    '''
    Return rotation matrices with rotational angle theta along an arbitary rotation axis.

    Parameters
    ----------
    theta : array-like
        Rotation angles of the matrix in radians.
    axis : ndarray, shape (N,3)
        Rotational axis for the rotation being performed, align with first index (i.e. axis[3] is the rotational axis
        for angle theta[3]).

    Returns
    -------
    rot_matrix : ndarray
        Corresponding rotation martrices for each input angle, align with the first index.

    See Also
    --------
    Rz : Return z-rotation matrices with rotational angle theta.
    Ry : Return y-rotation matrices with rotational angle theta.
    Rx : Return x-rotation matrices with rotational angle theta.
    '''
    if isinstance(theta, np.ndarray): # array version
        sin_theta = np.sin( theta )
        cos_theta = np.cos( theta )
        axis = axis / np.sqrt( np.sum(axis * axis, axis = 1) )[:,np.newaxis]
        rot_matrix = np.zeros( (theta.shape[0], 3, 3) )
        rot_matrix[:,0,0] = cos_theta + axis[:,0] ** 2 * (1 - cos_theta)
        rot_matrix[:,1,1] = cos_theta + axis[:,1] ** 2 * (1 - cos_theta)
        rot_matrix[:,2,2] = cos_theta + axis[:,2] ** 2 * (1 - cos_theta)
        rot_matrix[:,0,1] = axis[:,0] * axis[:,1] * (1 - cos_theta) - axis[:,2] * sin_theta
        rot_matrix[:,1,0] = axis[:,0] * axis[:,1] * (1 - cos_theta) + axis[:,2] * sin_theta
        rot_matrix[:,0,2] = axis[:,0] * axis[:,2] * (1 - cos_theta) + axis[:,1] * sin_theta
        rot_matrix[:,2,0] = axis[:,0] * axis[:,2] * (1 - cos_theta) - axis[:,1] * sin_theta
        rot_matrix[:,1,2] = axis[:,1] * axis[:,2] * (1 - cos_theta) - axis[:,0] * sin_theta
        rot_matrix[:,2,1] = axis[:,1] * axis[:,2] * (1 - cos_theta) + axis[:,0] * sin_theta
    else: #scalar
        sin_theta = math.sin( theta )
        cos_theta = math.cos( theta )
        axis = axis / math.sqrt( np.dot(axis, axis) )
        rot_matrix = np.identity(3)
        rot_matrix[0,0] = cos_theta + axis[0] ** 2 * (1 - cos_theta)
        rot_matrix[1,1] = cos_theta + axis[1] ** 2 * (1 - cos_theta)
        rot_matrix[2,2] = cos_theta + axis[2] ** 2 * (1 - cos_theta)
        rot_matrix[0,1] = axis[0] * axis[1] * (1 - cos_theta) - axis[2] * sin_theta
        rot_matrix[1,0] = axis[0] * axis[1] * (1 - cos_theta) + axis[2] * sin_theta
        rot_matrix[0,2] = axis[0] * axis[2] * (1 - cos_theta) + axis[1] * sin_theta
        rot_matrix[2,0] = axis[0] * axis[2] * (1 - cos_theta) - axis[1] * sin_theta
        rot_matrix[1,2] = axis[1] * axis[2] * (1 - cos_theta) - axis[0] * sin_theta
        rot_matrix[2,1] = axis[1] * axis[2] * (1 - cos_theta) + axis[0] * sin_theta
    return rot_matrix

def unitarize( R ):
    '''
    Enforce unitarity of the input matrix using Gram-Schmidt process.
    This function overwrites the input matrix.

    Parameters
    ----------
    R : ndarray, shape (N,3,3)
        (List of) input matrix to be unitarized.

    Returns
    -------
    R : ndarray, shape (N,3,3)
        Unitarized input matrix.
    '''
    if len( R.shape ) == 2: #1D case
        R[0] /= math.sqrt( R[0].dot( R[0]) )
        R[1] -= R[1].dot( R[0]) * R[0]
        R[1] /= math.sqrt(  R[1].dot( R[1]) )
        R[2] -= R[2].dot( R[0]) * R[0]
        R[2] -= R[2].dot( R[1]) * R[1]
        R[2] /= math.sqrt( R[2].dot( R[2]) )
    else: #2D case
        R[:,0] /= np.sqrt( np.sum( R[:,0] ** 2, axis=1) )[:,np.newaxis]
        R[:,1] -= np.einsum( 'ij,ij->i', R[:,1], R[:,0] )[:,np.newaxis] * R[:,0]
        R[:,1] /= np.sqrt( np.sum( R[:,1] ** 2, axis=1) )[:,np.newaxis]
        R[:,2] -= np.einsum( 'ij,ij->i', R[:,2], R[:,0] )[:,np.newaxis] * R[:,0]
        R[:,2] -= np.einsum( 'ij,ij->i', R[:,2], R[:,1] )[:,np.newaxis] * R[:,1]
        R[:,2] /= np.sqrt( np.sum( R[:,2] ** 2, axis=1) )[:,np.newaxis]
    return R

def params2coords( params ):
    '''
    Convert base-pair step parameters to bp-center coordinates and rotation matrix.

    Parameters
    ----------
    params : ndarray, shape (N,6)
        Input base-pair step parameters. Distances in unit of Angstroms, angles in unit of radians.
        Order = [Shift, Slide, Rise, Tilt, Roll, Twist]

    Returns
    -------
    o2 : ndarray, shape (N,3)
        Coordiantes for the bp-centers of the 2nd base-pair.
    R2 : ndarray, shape (N,3,3)
        Rotational matrix for the transformation.

    See Also
    --------
    coords2params : Convert base-pair-center coordinates and rotation matrix to base-pair step parameters.
    '''
    if len( params.shape ) == 1: #1D case
        Tau = math.sqrt( params[3] ** 2 + params[4] ** 2)
        Phi = math.atan2(params[3], params[4])
        o2 = Rz( params[5] * 0.5 - Phi ).dot( Ry( Tau * 0.5 ) ).dot( Rz( Phi ) ).dot( params[0:3] )
        R2 = Rz( -Phi + params[5] * 0.5 ).dot( Ry( Tau ) ).dot( Rz( params[5] * 0.5 + Phi ) )
    else: #2D case
        Tau = np.sqrt( params[:,3] ** 2 + params[:,4] ** 2)
        Phi = np.arctan2(params[:,3], params[:,4])
        o2 = np.einsum('ijk,ikl,ilm,im->ij', Rz( params[:,5] * 0.5 - Phi ), Ry( Tau * 0.5 ), Rz( Phi ), params[:,0:3] )
        R2 = np.einsum('ijk,ikl,ilm->ijm', Rz( -Phi + params[:,5] * 0.5 ), Ry( Tau ), Rz( params[:,5] * 0.5 + Phi ) )
    return o2, R2

def coords2params( o2, R2 ):
    '''
    Convert base-pair-center coordinates and rotation matrix to base-pair step parameters.
    The frame of bp1 is used as reference ( o1 = [0 0 0] and R1 = np.eye(3) )

    Parameters
    ----------
    o2 : ndarray, shape (N,3)
        Coordiantes for the bp-centers of the 2nd base-pair.
    R2 : ndarray, shape (N,3,3)
        Rotational matrix for the transformation.

    Returns
    -------
    params : ndarray, shape (N,6)
        Input base-pair step parameters. Distances in unit of Angstroms, angles in unit of radians.
        Order = [Shift, Slide, Rise, Tilt, Roll, Twist]

    See Also
    --------
    params2coords : Convert base-pair step parameters to bp-center coordinates and rotation matrix.
    '''
    if  len( o2.shape ) == 1: #1D case
        Tau = math.acos(R2[2,2])
        if math.sin( Tau ) > 0:
            alpha = math.atan2( R2[1,2],  R2[0,2])
            gamma = math.atan2( R2[2,1], -R2[2,0])
        else:
            alpha = math.atan2( -R2[1,2], -R2[0,2])
            gamma = math.atan2( -R2[2,1],  R2[2,0])
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
        d_vec = Rz( -Phi ).dot( Ry( -Tau * 0.5 ) ).dot( Rz( -Twst * 0.5 + Phi ) ).dot( o2 )
        return np.hstack( (d_vec, Tilt, Roll, Twst) )
    else: #2D case
        Tau = np.arccos(R2[:,2,2])
        sign_sin_Tau = np.sign( np.sin( Tau ) )
        alpha = np.arctan2( sign_sin_Tau * R2[:,1,2],  sign_sin_Tau * R2[:,0,2])
        gamma = np.arctan2( sign_sin_Tau * R2[:,2,1], -sign_sin_Tau * R2[:,2,0])
        Twst = alpha + gamma
        Phi = 0.5 * (gamma - alpha)
        Phi[Twst > pi] += pi
        Phi[Twst <= -pi] += pi
        Twst[Twst > pi] -= 2 * pi
        Twst[Twst <= -pi] += 2 *pi
        Roll = Tau * np.cos(Phi)
        Tilt = Tau * np.sin(Phi)
        d_vec = np.einsum('ijk,ikl,ilm,im->ij', Rz( -Phi ), Ry( -Tau * 0.5 ), Rz( -Twst * 0.5 + Phi ), o2 )
        return np.column_stack( (d_vec, Tilt, Roll, Twst) )

def params2data( params, frame0 = np.eye(3) ) :
    '''
    Convert list of base-pair step parameters to delta-r vectors and frames of each site.

    Parameters
    ----------
    params : ndarray, shape (N,6)
        Input base-pair step parameters.
        Order = [Shift, Slide, Rise, Tilt, Roll, Twist]
        Distance in unit of Angstroms, angle in unit of radians.
    frame0 : ndarray, shape (3,3), optional
        The frame for the 1st base-pair. Default set to np.eye(3).

    Returns
    -------
    dr : ndarray, shape (N,3)
        Delta-r vectors (dr[i] = r[i+1] - r[i]).
    frames : ndarray, shape (N+1,3,3)
        Frame of each base-pair.
    '''
    n_bp = params.shape[0] + 1
    frames = np.zeros( (n_bp,3,3) )
    o, R = params2coords( params )
    unitarize( R )
    for i in xrange(n_bp) :
        if i == 0 :
            frames[i] = unitarize( frame0 )
        else :
            frames[i] = unitarize( frames[i-1].dot(R[i-1]) )
    dr = np.einsum('ijk,ik->ij', frames[:i], o)
    return dr, frames

def dr2coord( dr ):
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
    '''
    coord = np.vstack( ( np.zeros(3), np.cumsum( dr, axis=0 ) ) )
    return coord

def ribbon_twist( dr, rb_vec, return_val_only=True ):
    '''
    Compute the ribbon-twist (i.e. supercoiling twist) of a helix.

    Parameters
    ----------
    dr : ndarray, shape (N,3)
        Input delta_r vectors (dr[i] = r[i+1] - r[i]).
    rb_vec : ndarray, shape (N+1,3)
        The ribbon vectors for the helix, equals to frame[:,:,1] in the current setting.
    return_val_only : bool, optional
        If True, return one float value for the overall twist of the system. Otherwise return the individual twists for each
        bp-step. Default set to True.

    Returns
    -------
    twist : float or ndarray of shape (N)
        The overall twist or the individual twists of the helix. See the 'return_val_only' parameter above.

    Raises
    ------
    ValueError
        If dr.shape[0] != rb_vec.shape[0] + 1

    See Also
    --------
    writhe_fuller : Compute the writhe of the helix using the Fuller's approximated single integral.
    writhe_exact : Compute the writhe of the helix using the exact Gauss double integral.
    '''
    if dr.shape[0] == rb_vec.shape[0] - 1:
        dr = np.vstack( (ez, dr, ez) )
    if dr.shape[0] != rb_vec.shape[0] + 1:
        raise ValueError('dr.shape[0] != rb_vec.shape[0] + 1')

    twist = np.empty( rb_vec.shape[0] - 1 )
    n = dr.shape[0]
    code = '''
    double norm, cos_angle, angle, sign;
    double v1[n-1][3];
    double alpha[n-1];
    for (int i = 0; i < (n-1); i++) {
        v1[i][0] = dr(i,1) * dr(i+1,2) - dr(i,2) * dr(i+1,1);
        v1[i][1] = dr(i,2) * dr(i+1,0) - dr(i,0) * dr(i+1,2);
        v1[i][2] = dr(i,0) * dr(i+1,1) - dr(i,1) * dr(i+1,0);
        if ( v1[i][0] == 0 && v1[i][1] == 0 && v1[i][2] == 0 ) {
            if (dr(i,2) != 0) {
                v1[i][0] = 1;
                v1[i][1] = 1;
                v1[i][2] = -dr(i,0) - dr(i,1);
            } else if (dr(i,1) != 0)  {
                v1[i][0] = 1;
                v1[i][2] = 1;
                v1[i][1] = -dr(i,0) - dr(i,2);
            } else {
                v1[i][1] = 1;
                v1[i][2] = 1;
                v1[i][0] = -dr(i,1) - dr(i,2);
            }
        }
        norm = sqrt( v1[i][0] * v1[i][0] + v1[i][1] * v1[i][1] + v1[i][2] * v1[i][2] );
        v1[i][0] /= norm;
        v1[i][1] /= norm;
        v1[i][2] /= norm;
    }

    for (int i = 0; i < (n-2); i++) {
        cos_angle = v1[i][0] * v1[i+1][0] + v1[i][1] * v1[i+1][1] + v1[i][2] * v1[i+1][2];
        angle;
        if (cos_angle > 1) {
            angle = pi;
        } else if (cos_angle < -1) {
            angle = -pi;
        } else {
            angle = acos( cos_angle );
        }
        sign = 0;
        sign += ( v1[i][1] * v1[i+1][2] - v1[i][2] * v1[i+1][1] ) * dr(i+1,0);
        sign += ( v1[i][2] * v1[i+1][0] - v1[i][0] * v1[i+1][2] ) * dr(i+1,1);
        sign += ( v1[i][0] * v1[i+1][1] - v1[i][1] * v1[i+1][0] ) * dr(i+1,2);
        if (sign > 0) {
            twist(i) = angle;
        } else {
            twist(i) = -angle;
        }
    }

    for (int i = 0; i < (n-1); i++) {
        cos_angle = v1[i][0] * rb_vec(i,0) + v1[i][1] * rb_vec(i,1) + v1[i][2] * rb_vec(i,2);
        angle;
        if (cos_angle > 1) {
            angle = pi;
        } else if (cos_angle < -1) {
            angle = -pi;
        } else {
            angle = acos( cos_angle );
        }
        sign = 0;
        sign += ( v1[i][1] * rb_vec(i,2) - v1[i][2] * rb_vec(i,1) ) * dr(i,0);
        sign += ( v1[i][2] * rb_vec(i,0) - v1[i][0] * rb_vec(i,2) ) * dr(i,1);
        sign += ( v1[i][0] * rb_vec(i,1) - v1[i][1] * rb_vec(i,0) ) * dr(i,2);
        if (sign > 0) {
            alpha[i] = angle;
        } else {
            alpha[i] = -angle;
        }
    }

    for (int i = 0; i < (n-2); i++) {
        twist(i) += alpha[i+1] - alpha[i];
        if (twist(i) > pi) {
            twist(i) -= 2 * pi;
        } else if (twist(i) < -pi) {
            twist(i) += 2 * pi;
        }
    }
    '''
    weave.inline( code, ['dr','rb_vec','twist','n','pi'], type_converters= weave.converters.blitz )
    if return_val_only :
        return np.sum( twist )
    else :
        return twist

def writhe_exact( dr ) :
    '''
    Compute the writhe of the helix using the exact Gauss double integral.
    O(N^2) complexity.

    Parameters
    ----------
    dr : ndarray, shape (N,3)
        Input delta_r vectors (dr[i] = r[i+1] - r[i]).

    Returns
    -------
    writhe : float
        The overall writhe of the helix.

    See Also
    --------
    writhe_fuller : Compute the writhe of the helix using the Fuller's approximated single integral.
    ribbon_twist : Compute the ribbon-twist of a helix.

    References
    ----------
    [1] Rossetto V, Maggs AC (2003) Writhing geometry of open DNA. J. Chem. Phys. 118: 9864-9874.
    [2] White JH (1969) Self-linking and the Gauss integral in higher dimensions. Am. J. Math. 91: 693-728.
    '''
    r0 = dr2coord( dr )
    r1 = r0[1:] - r0[0]
    r2 = r0[:-1] - r0[-1]

    r_size = r0.shape[0] - 1

    code = '''
    double writhe_int = 0;
    double norm, area, sign;
    double a[3];
    double v[4][3];
    for (int i = 0; i < r_size-2; ++i) {
        for (int j = i; j < r_size-2; ++j) {
            for (int k = 0; k < 3; ++k) {
                v[0][k] = r0(i,k) - r0(j+2,k);
                v[1][k] = r0(i,k) - r0(j+3,k);
                v[2][k] = r0(i+1,k) - r0(j+2,k);
                v[3][k] = r0(i+1,k) - r0(j+3,k);
            }
            for (int k = 0; k < 4; ++k) {
                norm = sqrt( v[k][0] * v[k][0] + v[k][1] * v[k][1] + v[k][2] * v[k][2] );
                v[k][0] /= norm;
                v[k][1] /= norm;
                v[k][2] /= norm;
            }

            //Area v0-v1-v3
            a[0] = v[0][0] * v[1][0] + v[0][1] * v[1][1] + v[0][2] * v[1][2];
            a[1] = v[0][0] * v[3][0] + v[0][1] * v[3][1] + v[0][2] * v[3][2];
            a[2] = v[1][0] * v[3][0] + v[1][1] * v[3][1] + v[1][2] * v[3][2];
            for (int k = 0; k < 3; ++k) {
                if (a[k] > 1) {
                    a[k] = 0.25 * pi;
                } else if (a[k] < -1) {
                    a[k] = -0.25 * pi;
                } else {
                    a[k] = 0.25 * acos( a[k] );
                }
            }
            area = tan( a[0] + a[1] + a[2] );
            area *= tan( a[2] - a[0] + a[1] );
            area *= tan( a[2] + a[0] - a[1] );
            area *= tan( a[0] + a[1] - a[2] );

            if (area < 0) {
                area = 0;
            } else {
                area = atan( sqrt( area ) ) * 4;
            }
            sign = 0;
            sign += ( v[3][1] * v[0][2] - v[3][2] * v[0][1] ) * v[1][0];
            sign += ( v[3][2] * v[0][0] - v[3][0] * v[0][2] ) * v[1][1];
            sign += ( v[3][0] * v[0][1] - v[3][1] * v[0][0] ) * v[1][2];


            if (sign > 0) {
                writhe_int += area;
            } else {
                writhe_int -= area;
            }

            //Area v0-v3-v2
            a[0] = v[0][0] * v[3][0] + v[0][1] * v[3][1] + v[0][2] * v[3][2];
            a[1] = v[0][0] * v[2][0] + v[0][1] * v[2][1] + v[0][2] * v[2][2];
            a[2] = v[3][0] * v[2][0] + v[3][1] * v[2][1] + v[3][2] * v[2][2];
            for (int k = 0; k < 3; ++k) {
                if (a[k] > 1) {
                    a[k] = 0.25 * pi;
                } else if (a[k] < -1) {
                    a[k] = -0.25 * pi;
                } else {
                    a[k] = 0.25 * acos( a[k] );
                }
            }
            area = tan( a[0] + a[1] + a[2] );
            area *= tan( a[2] - a[0] + a[1] );
            area *= tan( a[2] + a[0] - a[1] );
            area *= tan( a[0] + a[1] - a[2] );

            if (area < 0) {
                area = 0;
            } else {
                area = atan( sqrt( area ) ) * 4;
            }
            sign = 0;
            sign += ( v[0][1] * v[3][2] - v[0][2] * v[3][1] ) * v[2][0];
            sign += ( v[0][2] * v[3][0] - v[0][0] * v[3][2] ) * v[2][1];
            sign += ( v[0][0] * v[3][1] - v[0][1] * v[3][0] ) * v[2][2];

            if (sign > 0) {
                writhe_int += area;
            } else {
                writhe_int -= area;
            }
        }
    }
    return_val = writhe_int;
    '''
    writhe_internal = weave.inline( code, ['r_size','r0','pi'], type_converters= weave.converters.blitz )
    return writhe_internal + writhe_fuller( r1 ) + writhe_fuller( -r2 )


def writhe_fuller( dr, return_val_only=True ) :
    '''
    Compute the writhe of the helix using the Fuller's approximated single integral.
    Only guarantee to be correct modulo 2 pi. O(N) complexity.

    Parameters
    ----------
    dr : ndarray, shape = (N,3)
        Input delta_r vectors (dr[i] = r[i+1] - r[i]).
    return_val_only : bool, optional
        If True, return one float value for the overall writhe of the system. Otherwise return the individual writhe contribution
        for each bp-step. Default set to True.

    Returns
    -------
    writhe : float or ndarray of shape(N)
        The overall writhe of the helix or the individual writhe contribution for each bp-step. See the 'return_val_only' parameter above.

    See Also
    --------
    writhe_exact : Compute the writhe of the helix using the exact Gauss double integral.
    ribbon_twist : Compute the ribbon-twist of a helix.

    References
    ----------
    [1] Rossetto V, Maggs AC (2003) Writhing geometry of open DNA. J. Chem. Phys. 118: 9864-9874.
    [2] Fuller FB (1971) The writhing number of a space curve. PNAS 68: 815-819.
    '''

    n = dr.shape[0] - 1
    area = np.empty( n )
    code = '''
    double norm1, norm2, sign;
    double a[3];
    norm1 = sqrt( dr(0,0) * dr(0,0) + dr(0,1) * dr(0,1) + dr(0,2) * dr(0,2) );
    for (int i = 0; i < n; ++i) {
        norm2 = sqrt( dr(i+1,0) * dr(i+1,0) + dr(i+1,1) * dr(i+1,1) + dr(i+1,2) * dr(i+1,2) );
        a[0] = dr(i,2) / norm1;
        a[1] = dr(i+1,2) / norm2;
        a[2] = ( dr(i,0) * dr(i+1,0) + dr(i,1) * dr(i+1,1) + dr(i,2) * dr(i+1,2) ) / norm1 / norm2;
        norm1 = norm2;
        for (int j = 0; j < 3; ++j) {
            if (a[j] > 1) {
                a[j] = 0.25 * pi;
            } else if (a[j] < -1) {
                a[j] = -0.25 * pi;
            } else {
                a[j] = 0.25 * acos( a[j] );
            }
        }
        area(i) = tan( a[0] + a[1] + a[2] );
        area(i) *= tan( a[2] - a[0] + a[1] );
        area(i) *= tan( a[2] + a[0] - a[1] );
        area(i) *= tan( a[0] + a[1] - a[2] );

        if (area(i) < 0) {
            area(i) = 0;
        } else {
            area(i) = atan( sqrt( area(i) ) ) * 4;
        }
        sign = -dr(i,1) * dr(i+1,0) + dr(i,0) * dr(i+1,1);
        if (sign < 0) {
            area(i) = -area(i);
        }
    }
    '''
    weave.inline( code, ['dr','area','n','pi'], type_converters= weave.converters.blitz )

    if return_val_only :
        return np.sum(area)
    else :
        return area

def MC_acpt_rej( score_old, score_new, kT = kBT ):
    '''
    Decide whether to accept a Monte-Carlo move, given the old score and new score.

    Parameters
    ----------
    score_old : float
        score of the system before the proposed update.
    score_new : float
        score of the system after the proposed update.
    kT : float
        Temperature times Boltzmann constant.

    Returns
    -------
    out : bool
        True for accept, False for reject.
    '''
    exponent = (score_new - score_old) / kT
    if  exponent <= 1e-9 or math.exp(-exponent) >= random.random_sample():
        return True
    else:
        return False

def read_seq_from_fasta( fasta ):
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
        raise ValueError('The fasta file does not exist.')
    comment_symbols = ['>', ';', '#']
    seq = ''
    for line in open(fasta):
        if line[0] not in comment_symbols:
            seq += line.strip()
    return seq

#######################################
#Obsolete old functions
#def cos_regularize( v ) :
#    "Ensure the input values within the range of cos / sin ( -1 ~ 1 )"
#    v[ v > 1 ] = 1
#    v[ v < -1 ] = -1
#    return v
#
#def ribbon_twist_numpy( dr, rb_vec, return_val_only=True ):
#    if dr.shape[0] == rb_vec.shape[0] - 1:
#        dr = np.vstack( (ez, dr, ez) )
#    assert (dr.shape[0] == rb_vec.shape[0] + 1)
#
#    v1 = np.cross( dr[:-1], dr[1:] )
#    zero_entry = np.all( v1 == np.zeros(3), axis=1 )
#    if np.any(zero_entry):
#        dr1 = dr[:-1][zero_entry]
#        dr2 = np.random.rand( *dr1.shape )
#        v1[zero_entry] = np.cross( dr1, dr2 )
#
#    v1 /= np.sqrt( np.sum( v1 ** 2, axis=1) ) [:, np.newaxis]
#    angle1 = np.arccos( cos_regularize( np.einsum('ij,ij->i', v1[:-1], v1[1:]) ) )
#    angle1 *= np.sign( np.einsum('ij,ij->i', np.cross( v1[:-1], v1[1:] ), dr[1:-1] ) )
#
#    v2 = dr[:-1]
#    #v2 = (dr[:-1] + dr[1:]) * 0.5
#    #v2 /= np.sqrt( np.sum( v2 ** 2, axis=1) ) [:, np.newaxis]
#    #rb_vec = rb_vec - np.einsum('ij,ij->i', rb_vec, v2)[:,np.newaxis] * v2
#    #rb_vec = rb_vec / np.sqrt( np.sum( rb_vec ** 2, axis=1) ) [:, np.newaxis]
#
#    angle2 = np.arccos( cos_regularize( np.einsum('ij,ij->i', rb_vec, v1) ) )
#    angle2 *= np.sign( np.einsum('ij,ij->i', np.cross( v1, rb_vec ), v2 ) )
#    local_twist = np.mod( angle2[1:] - angle2[:-1] + angle1, 2 * pi )
#    local_twist[ local_twist > pi  ] -= 2 * pi
#
#    if return_val_only :
#        return np.sum( local_twist )
#    else :
#        return local_twist
#def writhe_fuller_numpy( v, return_val_only=True ) :
#    "Writhe calculation using Fuller's approximation. Correct to Mod 2. O(N)."
#
#    v_norm = v / np.sqrt( np.sum(  v ** 2, axis=1 ) )[:,np.newaxis]
#    v1 = v_norm[:-1]
#    v2 = v_norm[1:]
#
#    a = np.arccos( cos_regularize( v1[:,2] ) )#v0.dot(v1)
#    b = np.arccos( cos_regularize( v2[:,2] ) )#v0.dot(v2)
#    c = np.arccos( cos_regularize( np.sum( v1 * v2, axis=1 ) ) )
#    area = np.ones_like(a)
#    area *= np.tan( (a + b + c) * 0.25 )
#    area *= np.tan( (c - a + b) * 0.25 )
#    area *= np.tan( (c + a - b) * 0.25 )
#    area *= np.tan( (a + b - c) * 0.25 )
#    area[ area < 0 ] = 0
#    area = np.arctan( np.sqrt(area) ) * 4
#    v0_cross_v1 = np.column_stack( ( -v1[:,1], v1[:,0] ) )
#    area *= np.sign( np.sum( v0_cross_v1 * v2[:,:2], axis=1 ) )
#
#    if return_val_only :
#        return np.sum(area)
#    else :
#        return area
#
#def writhe_exact_numpy( v ) :
#    "The exact double integral formulation for writhe. O(N^2)."
#    r0 = dr2coord( v )
#    r1 = r0[1:] - r0[0]
#    r2 = r0[:-1] - r0[-1]
#
#    r_size = r0.shape[0]
#    r_mat = np.tile( r0[:,np.newaxis,:], (1,r_size,1) )
#    r_mat = r_mat - np.transpose( r_mat, (1,0,2) )
#    r_mat[ np.triu_indices(r_size,1) ] /= np.sqrt( np.sum( r_mat[ np.triu_indices(r_size,1) ] ** 2, axis=1 ) )[:,np.newaxis]
#
#    v0 = r_mat[:-1,:-1][ np.triu_indices(r_size-1,2) ]
#    v1 = r_mat[:-1,1:] [ np.triu_indices(r_size-1,2) ]
#    v2 = r_mat[1:,:-1] [ np.triu_indices(r_size-1,2) ]
#    v3 = r_mat[1:,1:]  [ np.triu_indices(r_size-1,2) ]
#
#    #Area v0-v1-v3
#    a = np.arccos( cos_regularize( np.einsum('ik,ik->i', v0, v1) ) )
#    b = np.arccos( cos_regularize( np.einsum('ik,ik->i', v0, v3) ) )
#    c = np.arccos( cos_regularize( np.einsum('ik,ik->i', v1, v3) ) )
#    new_area = np.ones_like(a)
#    new_area *= np.tan( (a + b + c) * 0.25 )
#    new_area *= np.tan( (c - a + b) * 0.25 )
#    new_area *= np.tan( (c + a - b) * 0.25 )
#    new_area *= np.tan( (a + b - c) * 0.25 )
#    new_area[ new_area < 0 ] = 0
#    new_area = np.arctan( np.sqrt(new_area) ) * 4
#    new_area *= np.sign( np.einsum('ik,ik->i', np.cross(v3,v0), v1) )
#    writhe_int = np.sum(new_area)
#
#    #Area v0-v3-v2
#    a = np.arccos( cos_regularize( np.einsum('ik,ik->i', v0, v3) ) )
#    b = np.arccos( cos_regularize( np.einsum('ik,ik->i', v0, v2) ) )
#    c = np.arccos( cos_regularize( np.einsum('ik,ik->i', v3, v2) ) )
#    new_area = np.ones_like(a)
#    new_area *= np.tan( (a + b + c) * 0.25 )
#    new_area *= np.tan( (c - a + b) * 0.25 )
#    new_area *= np.tan( (c + a - b) * 0.25 )
#    new_area *= np.tan( (a + b - c) * 0.25 )
#    new_area[ new_area < 0 ] = 0
#    new_area = np.arctan( np.sqrt(new_area) ) * 4
#    new_area *= np.sign( np.einsum('ik,ik->i', np.cross(v0,v3), v2) )
#    writhe_int += np.sum(new_area)
#
#    return writhe_int + writhe_fuller( r1 ) + writhe_fuller( -r2 )
###########################################

