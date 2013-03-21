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
from copy import copy
from __init__ import kBT

def wlc_bouchiat(Lp, L0, z, kT=kBT):
    r'''
    Worm-like chain fitting formula described in Bouchiat et al. 1999 paper.
    For fitting force-extension curve.

    Parameters
    ----------
    Lp : float
        Bending persistence length, in Angstrom.
    L0 : float
        Contour length, in Angstrom.
    z : float or 1D ndarray
        Average Z-extension of the helix, in Angstrom.
    kT : float, optional
        Boltzmann factor times temperature, default is at 298.15 K.

    Returns
    -------
    F : float or 1D ndarray
        Z-direction stretching force, in pN.

    See Also
    --------
    wlc_bouchiat_impl, wlc_bouchiat_impl_log : Implict worm-like chain fitting formula described in Bouchiat et al.
    wlc_bouchiat_impl_log : wlc_bouchiat_impl compared in log space.
    f_wlc_bouchiat_impl : Approximately solve the force from implicit wlc model by grid search.

    Notes
    -----
    Here is the function :
    .. math ::
       F=\frac{kT}{{{L}_{p}}}\left[ \frac{1}{4{{(1-z/{{L}_{0}})}^{2}}}-\frac{1}{4}+\frac{z}{{{L}_{0}}}+
       \sum\limits_{i=2}^{i\le 7}{{{\alpha }_{i}}{{\left( \frac{z}{{{L}_{0}}} \right)}^{i}}} \right]
    Where :math:`\alpha` = [-0.5164228, -2.737418, 16.07497, -38.87607, 39.49944, -14.17718] for i = 2~7.

    References
    ----------
    [1] Bouchiat C, Wang MD, Allemand J, Strick T, Block SM, et al. (1999) Estimating the persistence
    length of a worm-like chain molecule from force-extension measurements. Biophys. J. 76: 409-413.
    '''
    z_div_L = z / L0
    a = 0.25 / (1.0 - z_div_L) ** 2 - 0.25 + z_div_L
    alphas = [-0.5164228, -2.737418, 16.07497, -38.87607, 39.49944, -14.17718]
    b = copy(z_div_L)
    for alpha in alphas:
        b *= z_div_L
        a += alpha * b
    F = kT / Lp * a
    return F

def wlc_bouchiat_impl(Lp, L0, K0, z, F, kT=kBT):
    r'''
    Implicit worm-like chain fitting formula described in Bouchiat et al. 1999 paper.
    For fitting force-extension curve at high-force.

    Parameters
    ----------
    Lp : float
        Bending persistence length, in Angstrom.
    L0 : float
        Contour length, in Angstrom.
    K0 : float
        Stretch modulus, in pN.
    z : float or 1D ndarray
        Average Z-extension of the helix, in Angstrom.
    F : float or  1D ndarray
        Z-direction stretching force, in pN.
    kT : float, optional
        Boltzmann factor times temperature, default is at 298.15 K.

    Returns
    -------
    zero : float or 1D ndarray
        Difference value, should approach zero at perfect fit.

    See Also
    --------
    wlc_bouchiat : Worm-like chain fitting formula described in Bouchiat et al.
    wlc_bouchiat_impl_log : wlc_bouchiat_impl compared in log space.
    f_wlc_bouchiat_impl : Approximately solve the force from implicit wlc model by grid search.

    Notes
    -----
    Equation :
    .. math ::
       \frac{kT}{{{L}_{p}}}\left[ \frac{1}{4{{(1-l)}^{2}}}-\frac{1}{4}+l+
       \sum\limits_{i=2}^{i\le 7}{{{\alpha }_{i}}{{l}^{i}}} \right]-F=0
    With :math:`l=z/{{L}_{0}}-F/{{K}_{0}}`

    References
    ----------
    [1] Bouchiat C, Wang MD, Allemand J, Strick T, Block SM, et al. (1999) Estimating the persistence
    length of a worm-like chain molecule from force-extension measurements. Biophys. J. 76: 409-413.
    '''
    l = z / L0 - F / K0
    a = 0.25 / (1.0 - l) ** 2 - 0.25 + l
    alphas = [-0.5164228, -2.737418, 16.07497, -38.87607, 39.49944, -14.17718]
    b = copy(l)
    for alpha in alphas:
        b *= l
        a += alpha * b
    zero = kT / Lp * a - F
    return zero

def wlc_bouchiat_impl_log(Lp, L0, K0, z, F, kT=kBT):
    r'''
    Similar to `wlc_bouchiat_impl`, but comparison is done in log10(F) space instead.
    In practice it is more robust. See `wlc_bouchiat_impl` for details.

    See Also
    --------
    wlc_bouchiat : Worm-like chain fitting formula described in Bouchiat et al.
    wlc_bouchiat_impl : Implict worm-like chain fitting formula described in Bouchiat et al.
    f_wlc_bouchiat_impl : Approximately solve the force from implicit wlc model by grid search.
    '''
    l = z / L0 - F / K0
    a = 0.25 / (1.0 - l) ** 2 - 0.25 + l
    alphas = [-0.5164228, -2.737418, 16.07497, -38.87607, 39.49944, -14.17718]
    b = copy(l)
    for alpha in alphas:
        b *= l
        a += alpha * b
    zero = np.log10(kT / Lp * a) - np.log10(F)
    return zero

def f_wlc_bouchiat_impl(Lp, L0, K0, z, kT=kBT):
    '''
    Approximately solve the force from implicit wlc model by grid search.

    Parameters
    ----------
    Lp : float
        Bending persistence length, in Angstrom.
    L0 : float
        Contour length, in Angstrom.
    K0 : float
        Stretch modulus, in pN.
    z : float or 1D ndarray
        Average Z-extension of the helix, in Angstrom.
    kT : float, optional
        Boltzmann factor times temperature, default is at 298.15 K.

    Returns
    -------
    F : float or 1D ndarray
        Z-direction stretching force, in pN.

    See Also
    --------
    wlc_bouchiat : Worm-like chain fitting formula described in Bouchiat et al.
    wlc_bouchiat_impl : Implict worm-like chain fitting formula described in Bouchiat et al.
    wlc_bouchiat_impl_log : wlc_bouchiat_impl compared in log space.
    '''
    F_list= np.exp( np.linspace( np.log(0.001), np.log(100), 30000 ) )
    if isinstance(z,float): #float version
        zero = wlc_bouchiat_impl(Lp, L0, K0, z, F_list, kT=kBT)
        return F_list[np.argmin( np.abs(zero) )]
    else: #ndarray version
        F = []
        for z_indv in z:
            zero = wlc_bouchiat_impl(Lp, L0, K0, z_indv, F_list, kT=kBT)
            F.append( F_list[np.argmin( np.abs(zero) )] )
        return np.array(F)

def moroz_3rd(Lp, C, F, kT=kBT):
    r'''
    3rd order Moroz-Nelson function for effective torsional persistence.

    Parameters
    ----------
    Lp : float
        Bending persistence length, in Angstrom.
    C : float
        Torsional persistence length, in Angstrom.
    F : float or 1D ndarray
        Z-direction stretching force, in pN.

    Returns
    -------
    Ceff : float or 1D ndarray
        Effective torsional persistence length in Angstrom. Ceff = L0 / Var(Lk).

    See Also
    --------
    moroz_1st : 1st order Moroz-Nelson function for effective torsional persistence.

    Notes
    -----
    Equation :
    .. math::
       {{C}_{eff}}=C\left( 1-\frac{M}{4K}+\frac{{{M}^{2}}-2M}{16{{K}^{2}}}
       -\frac{4{{M}^{3}}-16{{M}^{2}}+21M}{256{{K}^{3}}} \right) \\
       M=\frac{C}{{{L}_{p}}},K=\sqrt{\frac{{{L}_{p}}F}{kT}}

    References
    ----------
    [1] Moroz JD, Nelson P (1997) Torsional directed walks, entropic
    elasticity, and DNA twist stiffness. PNAS 94: 14418-14422.
    [2] Moroz JD, Nelson P (1998) Entropic Elasticity of Twist-Storing
    Polymers. Macromolecules 31: 6333-6347.
    [3] Lipfert J, Wiggin M, Kerssemakers JWJ, Pedaci F, Dekker NH (2011)
    Freely orbiting magnetic tweezers to directly monitor changes in the twist
    of nucleic acids. Nat. Comm. 2: 439.
    '''
    K = np.sqrt( Lp * F / kT )
    M = C / Lp
    Ceff = ( 1 - M / K / 4 + ( M ** 2 - 2 * M ) / (K ** 2) / 16 -
        ( 4 * M ** 3 - 16 * M ** 2 + 21 * M ) / (K ** 3) / 256 ) * C
    return Ceff

def moroz_1st(Lp, C, F, kT=kBT):
    '''
    1st order Moroz-Nelson function for effective torsional persistence.
    See `moroz_3rd` for detailed explanation.

    See Also
    --------
    moroz_3rd : 3rd order Moroz-Nelson function for effective torsional persistence.
    '''
    K = np.sqrt( Lp * F / kT )
    M = C / Lp
    Ceff = ( 1 - M / K / 4 ) * C
    return Ceff
