#!/usr/bin/env python
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
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
from util import params2data, writhe_exact, writhe_fuller, ribbon_twist, params2coords, unitarize, dr2coord
from __init__ import ez

#####Pose object storing the state info of the helix#####
class HelixPose(object):
    '''
    Pose object storing the state info of the helix.

    Parameters
    ----------
    input_file : str, optional
        Input file (.npz) that stores the helix pose data. Overides 'init_params' if both specified.
    init_params : ndarray, shape (6), optional
        Parameter set used for initializing the helix pose. The initialized pose will have the same bp-step parameters throughout.
    n_bp : int, optional
        Number of base-pairs in the model. Must be specified if init_params is given.
    compute_tw_wr : bool, optional
        Whether to compute twist and writhe (Fuller's writhe) during system update. Should be set to True for
        link-constrained simulations.

    Attributes
    ----------
    compute_tw_wr : bool
        See Parameters section above.
    writhe_exact : float
        Exact writhe of the helix.
    writhe_fuller : float
        Fuller's writhe of the helix.
    twist : float
        Supercoiling twist of the helix.
    coord_terminal : ndarray
        Coordinate (x,y,z) of the center of last base-pair.
    frame_terminal : ndarray
        Coordinate frame of the last base-pair. Each column represent the cooresponding axis (frame[:,0] is the x-axis etc.)
    z_terminal : float
        Z-component of coord_terminal
    link_fuller : float
        Link of the helix computed using Fuller's approximation.
    link_exact : float
        Exact link of the helix.
    coord : ndarray
        Coordinates of all base-pairs in the helix.
    dr : ndarray
        Delta-r vectors of the entire helix.
    frame : ndarray
        Frames of all bp-steps in the helix.
    params : ndarray
        List of all bp-step parameters in the helix.
    n_bp : int
        Number of base-pairs in the helix.
    rb_vec : ndarray
        Ribbon vectors of all base-pairs in the helix.

    Raises
    ------
    ValueError
        If n_bp < 2
    '''
    def __init__(self, input_file=None, init_params=None, n_bp=0, compute_tw_wr=False):
        #core data
        self._n_bp = None
        self._dr = None
        self._frames = None
        self._params = None
        self._writhe_data = None
        self._twist_data = None
        self._compute_tw_wr = False

        #Auxillary cached data
        self._obs = {}
        self._obs['twist'] = None
        self._obs['writhe_fuller'] = None
        self._obs['writhe_exact'] = None
        self._obs['coord_terminal'] = None
        self._obs['frame_terminal'] = None
        self._obs['score'] = None

        if input_file != None:
            self.load_from_file(input_file)
            self.compute_tw_wr = compute_tw_wr
        elif init_params != None:
            if n_bp < 2:
                raise ValueError('n_bp < 2.')
            self._n_bp = n_bp
            self._params = np.tile( init_params, (n_bp-1,1) ).copy()
            self._dr, self._frames = params2data( self._params )
            self.compute_tw_wr = compute_tw_wr

    def load_from_file(self, input_file):
        '''
        Load pose data from an input file.

        Parameters
        ----------
        input_file : str
            Name of the input file containing the pose data.
        '''
        data = np.load(input_file)
        self._params = data['params']
        self._dr = data['dr']
        self._frames = data['frames']
        self._n_bp = self._frames.shape[0]
        if self._n_bp < 2:
            raise ValueError('n_bp < 2.')

    def copy(self):
        '''
        Return a copy of the current pose.

        Returns
        -------
        pose_copy : HelixPose
            Copy of the current pose.
        '''
        pose_copy = HelixPose()
        pose_copy._n_bp = self._n_bp
        pose_copy._dr = self._dr.copy()
        pose_copy._frames = self._frames.copy()
        pose_copy._params = self._params.copy()
        pose_copy._writhe_data = self._writhe_data.copy()
        pose_copy._twist_data = self._twist_data.copy()
        pose_copy._obs = self._obs.copy()
        pose_copy._compute_tw_wr = self._compute_tw_wr
        return pose_copy

    def write2disk(self, filename):
        '''
        Write the current pose to disk.

        Parameters
        ----------
        filename : str
            Name of the output file.
        '''
        np.savez( filename, params=self._params, dr=self._dr, frames=self._frames )

    @property
    def compute_tw_wr(self):
        return self._compute_tw_wr

    @compute_tw_wr.setter
    def compute_tw_wr(self, val):
        if val:
            self._compute_tw_wr = True
            self._update_all_writhe()
            self._update_all_twist()
        else:
            self._compute_tw_wr = False
    #######################
    #Update system
    def _update_all_writhe(self):
        "Update Writhe data with Fuller's approximation."
        self._writhe_data = writhe_fuller( self._dr, return_val_only=False )

    def _update_all_twist(self):
        "Update Twist data."
        self._twist_data = ribbon_twist( self._dr, self._frames[:,:,1], return_val_only=False )

    def _update_indv_twist(self, full_update=False):
        '''
        Update the twist for individual bp-step.

        Parameters
        ----------
        full_update : bool, optional
            True for full update and False for trial update.
        '''
        if full_update:
            twist = self.twist_data
            dr = self._dr
            self._mat = np.eye(3)
        else:
            twist = self._twist_data.copy()
            dr = self._dr_new
        i = self._working_bp_step
        #Special case
        if self._n_bp <= 4:
            frames = self._frames.copy()
            frames[(i+1):] = np.einsum('jk,ikl->ijl', self._mat, frames[(i+1):])
            twist[:] = ribbon_twist( self._dr, frames[:,:,1], False )
            return twist

        if i == 0 or i == 1:
            frames = self._frames[:(i+3)].copy()
            frames[(i+1):] = np.einsum('jk,ikl->ijl', self._mat, frames[(i+1):])
            twist[:(i+2)] = ribbon_twist( np.vstack( (ez, dr[:(i+3)]) ), frames[:,:,1], False )
            frames = np.einsum('jk,ikl->ijl', self._mat, self._frames[-2:])
            twist[-1] = ribbon_twist( np.vstack( (dr[-2:], ez) ), frames[:,:,1], False )
        elif i == self._n_bp - 2 or i == self._n_bp - 3:
            frames = self._frames[(i-1):].copy()
            frames[2:] = np.einsum('jk,ikl->ijl', self._mat, frames[2:])
            twist[(i-1):] = ribbon_twist( np.vstack( (dr[(i-2):], ez) ), frames[:,:,1], False )
        else:
            frames = self._frames[(i-1):(i+3)].copy()
            frames[2:] = np.einsum('jk,ikl->ijl', self._mat, frames[2:])
            twist[(i-1):(i+2)] = ribbon_twist( dr[(i-2):(i+3)], frames[:,:,1], False )
            frames = np.einsum('jk,ikl->ijl', self._mat, self._frames[-2:])
            twist[-1] = ribbon_twist( np.vstack( (dr[-2:], ez) ), frames[:,:,1], False )
        return twist

    def update(self, i, params, o=None, R=None ):
        '''
        Full update of the i-th bp-step.

        Parameters
        ----------
        i : int
            The index of the bp-step to be updated.
        params : ndarray
            Input base-pair step parameter set.
        o : ndarray, optional
            Coordiantes for the bp-centers of the 2nd base-pair.
        R : ndarray, optional
            Rotational matrix for the transformation.

        Raises
        ------
        ValueError
            If i < 0 or i >= n_bp.
        '''
        if i < 0 or i >= self._n_bp:
            raise ValueError('i < 0 or i >= n_bp.')
        if o == None or R == None:
            o, R = params2coords( params )
        self._obs_clear()

        mat = unitarize( self._frames[i].dot(R).dot(self._frames[i+1].T) )
        dr1 = self._frames[i].dot(o)
        self._frames[(i+1):] = np.einsum('jk,ikl->ijl', mat, self._frames[(i+1):])
        self._dr[i:] = np.vstack( ( dr1, np.einsum('jk,ik->ij', mat, self._dr[(i+1):]) ) )
        self._params[i] = params
        if self.compute_tw_wr:
            self._working_bp_step = i
            self._update_indv_twist(True)
            if i == 0 or i == 1:
                self._update_all_writhe()
            else:
                self._writhe_data[(i-1):] = writhe_fuller(self._dr[(i-1):], return_val_only=False)


    def update_trial(self, i, params, o=None, R=None ):
        '''
        Trail update of the i-th bp-step.
        Following accept_update or reject_update is required.

        Parameters
        ----------
        i : int
            The index of the bp-step to be updated.
        params : ndarray
            Input base-pair step parameter set.
        o : ndarray, optional
            Coordiantes for the bp-centers of the 2nd base-pair.
        R : ndarray, optional
            Rotational matrix for the transformation.

        Raises
        ------
        ValueError
            If i < 0 or i >= n_bp.
        '''
        if i < 0 or i >= self._n_bp:
            raise ValueError('i < 0 or i >= n_bp.')
        if o == None or R == None:
            o, R = params2coords( params )

        self._params_new = params
        self._mat = unitarize( self._frames[i].dot(R).dot(self._frames[i+1].T) )
        self._dr1 = self._frames[i].dot(o)
        self._working_bp_step = i
        self._obs_bkup()

        if self.compute_tw_wr:
            self._dr_new = np.vstack( ( self._dr[:i], self._dr1, np.einsum('jk,ik->ij', self._mat, self._dr[(i+1):]) ) )
            self._twist_new = self._update_indv_twist()
            if i == 0 or i == 1:
                self._writhe_new = writhe_fuller( self._dr_new, False )
            else:
                self._writhe_new = np.hstack((self._writhe_data[:(i-1)], writhe_fuller(self._dr_new[(i-1):], False)))
            self._obs['frame_terminal'] = self._frames[-1].dot( self._mat )
            self._obs['twist'] = np.sum(self._twist_new)
            self._obs['writhe_fuller'] = np.sum(self._writhe_new)
            self._obs['coord_terminal'] = np.sum(self._dr_new, axis=0)
            updated_keys = ['writhe_fuller', 'coord_terminal', 'twist', 'frame_terminal']
        else:
            self._obs['coord_terminal'] = np.sum(self._dr[:i], axis=0) + self._dr1 + self._mat.dot( np.sum(self._dr[(i+1):], axis=0) )
            self._obs['frame_terminal'] = self._frames[-1].dot( self._mat )
            updated_keys = ['coord_terminal', 'frame_terminal']
        self._obs_clear( exclude=updated_keys )

    def accept_update(self):
        'Accept a trial update.'
        mat = self._mat
        i = self._working_bp_step
        if self.compute_tw_wr:
            self._frames[(i+1):] = np.einsum('jk,ikl->ijl', mat, self._frames[(i+1):])
            self._dr = self._dr_new
            self._params[i] = self._params_new
            self._writhe_data = self._writhe_new
            self._twist_data = self._twist_new
        else:
            self._frames[(i+1):] = np.einsum('jk,ikl->ijl', mat, self._frames[(i+1):])
            self._dr[i:] = np.vstack( ( self._dr1, np.einsum('jk,ik->ij', mat, self._dr[(i+1):]) ) )
            self._params[i] = self._params_new

    def reject_update(self):
        'Reject a trial update.'
        self._obs_revert()

    ########################
    #Cleanup and backup observables
    def _obs_clear(self, exclude=None):
        ''''
        Clear self._obs.

        Parameters
        ----------
        exclude : list, optional
            Keys that are excluded from the clear action.
        '''
        if exclude == None:
            for j in self._obs.iterkeys():
                self._obs[j] = None
        else:
            for j in self._obs.iterkeys():
                if j not in exclude:
                    self._obs[j] = None

    def _obs_bkup(self):
        'Backup the current _obs.'
        self._obs_old = self._obs.copy()

    def _obs_revert(self):
        'Revert to the backup _obs.'
        self._obs = self._obs_old
    ########################
    #Returning observables
    @property
    def writhe_exact(self):
        if self._obs['writhe_exact'] == None:
             self._obs['writhe_exact']  = writhe_exact(self._dr)
        return self._obs['writhe_exact']

    @property
    def writhe_fuller(self):
        if self._obs['writhe_fuller'] == None:
            if not self.compute_tw_wr:
                self._obs['writhe_fuller'] = writhe_fuller( self._dr )
            else:
                self._obs['writhe_fuller'] = np.sum(self._writhe_data)
        return self._obs['writhe_fuller']

    @property
    def twist(self):
        if self._obs['twist'] == None:
            if not self.compute_tw_wr:
                self._obs['twist'] = ribbon_twist( self._dr, self._frames[:,:,1] )
            else:
                self._obs['twist'] = np.sum(self._twist_data)
        return self._obs['twist']

    @property
    def coord_terminal(self):
        if self._obs['coord_terminal'] == None:
             self._obs['coord_terminal'] = np.sum( self._dr, axis=0 )
        return self._obs['coord_terminal']

    @property
    def frame_terminal(self):
        if self._obs['frame_terminal'] == None:
             self._obs['frame_terminal'] = self._frames[-1].copy()
        return self._obs['frame_terminal']

    @property
    def z_terminal(self):
        return self.coord_terminal[2]

    @property
    def link_fuller(self):
        return self.twist + self.writhe_fuller

    @property
    def link_exact(self):
        return self.twist + self.writhe_exact

    @property
    def coord(self):
        return dr2coord( self._dr )

    @property
    def dr(self):
        return self._dr.copy()

    @property
    def frames(self):
        return self._frames.copy()

    @property
    def params(self):
        return self._params.copy()

    @property
    def n_bp(self):
        return self._n_bp

    @property
    def rb_vec(self):
        return self._frames[:,:,1].copy()
    #########################
    #Plotting functions
    def plot_helix(self, rb_width=16.0, color='kb', show=True, fig_ax=None) :
        '''
        Plot the helix using matplotlib + mplot3d.

        Parameters
        ----------
        rb_width : float, optional
            Width of the helix ribbon in Angstrom.
        color : str, optional
            Colors of the plot (matplotlib color codes), first letter for backbond and second letter for base-pair.
        show : bool, optional
            Whether to invoke pyplot.show() method at the end.
        fig_ax : mplot3d.Axes3D, optional
            Input mplot3d.Axes3D object. For ploting multiple helix on the same graph.
        '''
        rb_vec = self.rb_vec
        coord = self.coord
        bb1 = coord - rb_vec * rb_width * 0.5
        bb2 = coord + rb_vec * rb_width * 0.5
        if fig_ax == None:
            fig = plt.figure()
            ax = p3.Axes3D(fig)
        else:
            ax = fig_ax
        ax.plot3D(bb1[:,0], bb1[:,1], bb1[:,2], color[0]+'-')
        ax.plot3D(bb2[:,0], bb2[:,1], bb2[:,2], color[0]+'-')
        for i in xrange( bb1.shape[0] ):
            ax.plot3D( np.array( [bb1[i,0], bb2[i,0]] ),
                       np.array( [bb1[i,1], bb2[i,1]] ),
                       np.array( [bb1[i,2], bb2[i,2]] ), color[1]+'-')

        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        if show:
            plt.show()

    def plot_centerline(self, color='k', show=True, fig_ax=None) :
        '''
        Plot the helix center-line using matplotlib + mplot3d.

        Parameters
        ----------
        color : str, optional
            Colors of the centerline (matplotlib color codes).
        show : bool, optional
            Whether to invoke pyplot.show() method at the end.
        fig_ax : mplot3d.Axes3D, optional
            Input mplot3d.Axes3D object. For ploting multiple helix on the same graph.
        '''
        coord = self.coord
        if fig_ax == None:
            fig = plt.figure()
            ax = p3.Axes3D(fig)
        else:
            ax = fig_ax
        ax.plot3D(coord[:,0], coord[:,1], coord[:,2], color+'-')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        if show:
            plt.show()
#######################################
