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
    Assumes the chain has constant length (inextensive model).
    For fitting force-extension curve in low-to-medium force regime (< 10 pN).

    Parameters
    ----------
    Lp : float
        Bending persistence length, in Angstrom.
    L0 : float
        Contour length, in Angstrom.
    z : float or 1D ndarray
        Average Z-extension of the helix, in Angstrom.
    kT : float, optional
        Temperature times Boltzmann constant, in pN.A.

    Returns
    -------
    F : float or 1D ndarray
        Z-direction stretching force, in pN.

    See Also
    --------
    wlc_bouchiat_impl : Implict worm-like chain fitting formula described in Bouchiat et al.
    f_wlc_bouchiat_impl : Approximately solve the force from implicit wlc model by grid search.

    Notes
    -----
    The fitting function [BW]_ is

    .. math::
       F = \frac{kT}{L_p} \left[ \frac{1}{4 (1 - z / L_0)^2} - \frac{1}{4} + \frac{z}{L_0} +
       \sum\limits_{i=2}^{i \le 7}{ \alpha_i {\left( \frac{z}{L_0}  \right)}^i} \right]

    Where :math:`\alpha` = [-0.5164228, -2.737418, 16.07497, -38.87607, 39.49944, -14.17718] for i = 2~7.

    References
    ----------
    .. [BW] Bouchiat C, Wang MD, Allemand J, Strick T, Block SM, et al. (1999) Estimating the persistence
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

def wlc_bouchiat_impl(Lp, L0, S, z, F, kT=kBT):
    r'''
    Implicit worm-like chain fitting formula described in Bouchiat et al. 1999 paper.
    Comparison is done in log10(F) space instead (more robust in practice).
    An extensive rod factor is added for fitting force-extension curve at high-force.

    Parameters
    ----------
    Lp : float
        Bending persistence length, in Angstrom.
    L0 : float
        Contour length, in Angstrom.
    S : float
        Stretch modulus, in pN.
    z : float or 1D ndarray
        Average Z-extension of the helix, in Angstrom.
    F : float or  1D ndarray
        Z-direction stretching force, in pN.
    kT : float, optional
        Temperature times Boltzmann constant, in pN.A.

    Returns
    -------
    zero : float or 1D ndarray
        Difference value, should approach zero at perfect fit.

    See Also
    --------
    wlc_bouchiat : Worm-like chain fitting formula described in Bouchiat et al.
    f_wlc_bouchiat_impl : Approximately solve the force from implicit wlc model by grid search.

    Notes
    -----
    The fitting function [BW]_ is

    .. math::
       \log_{10} \left[  \frac{kT}{L_p} \left( \frac{1}{4(1-l)^2} - \frac{1}{4} + l +
       \sum \limits_{i=2}^{i \le 7}{\alpha_i l^i} \right) \right]- \log_{10} (F) = 0

    With :math:`l = z / L_0 - F / S`

    See `wlc_bouchiat` for reference papers.
    '''
    l = z / L0 - F / S
    a = 0.25 / (1.0 - l) ** 2 - 0.25 + l
    alphas = [-0.5164228, -2.737418, 16.07497, -38.87607, 39.49944, -14.17718]
    b = copy(l)
    for alpha in alphas:
        b *= l
        a += alpha * b
    zero = np.log10(kT / Lp * a) - np.log10(F)
    return zero

def f_wlc_bouchiat_impl(Lp, L0, S, z, kT=kBT):
    '''
    Approximately solve the force from implicit wlc model by grid search.

    Parameters
    ----------
    Lp : float
        Bending persistence length, in Angstrom.
    L0 : float
        Contour length, in Angstrom.
    S : float
        Stretch modulus, in pN.
    z : float or 1D ndarray
        Average Z-extension of the helix, in Angstrom.
    kT : float, optional
        Temperature times Boltzmann constant, in pN.A.

    Returns
    -------
    F : float or 1D ndarray
        Z-direction stretching force, in pN.

    See Also
    --------
    wlc_bouchiat : Worm-like chain fitting formula described in Bouchiat et al.
    wlc_bouchiat_impl : Implict worm-like chain fitting formula described in Bouchiat et al.
    '''
    F_list= np.exp( np.linspace( np.log(0.001), np.log(100), 30000 ) )
    if isinstance(z,float): #float version
        zero = wlc_bouchiat_impl(Lp, L0, S, z, F_list, kT=kBT)
        return F_list[np.argmin( np.abs(zero) )]
    else: #ndarray version
        F = []
        for z_indv in z:
            zero = wlc_bouchiat_impl(Lp, L0, S, z_indv, F_list, kT=kBT)
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
    The fitting function, orginally proposed in reference [M1]_ [M2]_,
    and summarized in [LW]_, is

    .. math::
       C_{eff} = C \left( 1-\frac{M}{4K} + \frac{M^2-2M}{16K^2}
       - \frac{4M^3-16M^2+21M}{256K^3} \right )

    With

    .. math::
       M = \frac{C}{L_p}, \quad K=\sqrt{ \frac{L_p F}{kT} }

    References
    ----------
    .. [M1] Moroz JD, Nelson P (1997) Torsional directed walks, entropic
       elasticity, and DNA twist stiffness. PNAS 94: 14418-14422.

    .. [M2] Moroz JD, Nelson P (1998) Entropic Elasticity of Twist-Storing
       Polymers. Macromolecules 31: 6333-6347.

    .. [LW] Lipfert J, Wiggin M, Kerssemakers JWJ, Pedaci F, Dekker NH (2011)
       Freely orbiting magnetic tweezers to directly monitor changes in the twist
       of nucleic acids. Nat. Comm. 2: 439.
    '''
    K = np.sqrt( Lp * F / kT )
    M = C / Lp
    Ceff = ( 1 - M / K / 4 + ( M ** 2 - 2 * M ) / (K ** 2) / 16 -
        ( 4 * M ** 3 - 16 * M ** 2 + 21 * M ) / (K ** 3) / 256 ) * C
    return Ceff

def moroz_1st(Lp, C, F, kT=kBT):
    r'''
    1st order Moroz-Nelson function for effective torsional persistence.
    See `moroz_3rd` for detailed explanation.

    See Also
    --------
    moroz_3rd : 3rd order Moroz-Nelson function for effective torsional persistence.

    Notes
    -----
    The fitting function is

    .. math::
       C_{eff} = C \left( 1-\frac{M}{4K} \right)

    With

    .. math::
       M = \frac{C}{L_p}, \quad K=\sqrt{ \frac{L_p F}{kT} }
    '''
    K = np.sqrt( Lp * F / kT )
    M = C / Lp
    Ceff = ( 1 - M / K / 4 ) * C
    return Ceff
