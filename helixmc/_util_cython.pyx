#cython: boundscheck=False
#cython: wraparound=False

#Note: here we used Typed Memoryviews in Cython

import numpy as np
cimport numpy as np
cimport cython
from cpython cimport bool
from __init__ import ez

cdef extern from "math.h" nogil:
    double sqrt(double x)
    double asin(double x)
    double acos(double x)

# Useful constants and vectors #
cdef double pi = np.pi

cdef double ex_c[3]
ex_c[0] = 1
ex_c[1] = 0
ex_c[2] = 0
#####################################################################
# Useful C utility functions #
cdef inline double arcsin(double x) nogil:
    "Arcsin with tolerence of x > 1 or x < -1 (due to numeric error)"
    global pi
    if x >= 1:
        return 0.5 * pi
    elif x <= -1:
        return -0.5 * pi
    else:
        return asin(x)

cdef inline double arccos(double x) nogil:
    "Arccos with tolerence of x > 1 or x < -1 (due to numeric error)"
    global pi
    if x >= 1:
        return 0
    elif x <= -1:
        return pi
    else:
        return acos(x)

cdef inline void cross_ez(double * r1, double * rc) nogil:
    "rc = cross(r1, ez)"
    rc[0] = r1[1]
    rc[1] = -r1[0]
    rc[2] = 0

cdef inline void cross_ez_rev(double * r1, double * rc) nogil:
    "rc = cross(ez, r1)"
    rc[0] = -r1[1]
    rc[1] = r1[0]
    rc[2] = 0

cdef inline void cross(double * r1, double * r2, double * rc) nogil:
    "rc = cross(r1, r2)"
    rc[0] = r1[1] * r2[2] - r1[2] * r2[1]
    rc[1] = r1[2] * r2[0] - r1[0] * r2[2]
    rc[2] = r1[0] * r2[1] - r1[1] * r2[0]

cdef inline double dot(double * r1, double * r2) nogil:
    "Dot product"
    return r1[0] * r2[0] + r1[1] * r2[1] + r1[2] * r2[2]


@cython.cdivision(True)
cdef inline void normalize(double * r) nogil:
    "Normalize the input vector"
    cdef double norm = sqrt(dot(r, r))
    r[0] /= norm
    r[1] /= norm
    r[2] /= norm

cdef inline bint is_vec_nil(double * r) nogil:
    "Check if the vector is all zero."
    if r[0] == 0 and r[1] == 0 and r[2] == 0:
        return True
    else:
        return False

cdef inline void normal_vec(double * r, double * rp) nogil:
    "Compute an arbitary normal vector"
    global ex_c
    cross(r, ex_c, rp)  # cross(r, ex)
    if is_vec_nil(rp):
        rp[0] = 0
        rp[1] = 0
        rp[2] = 1  # Since r is parallel to ex, just take ez as normal vector

##########################################################
# Compute the internal writhe in pure Cython #
cdef double _writhe_int(double[:, ::1] r) nogil:
    "Compute the internal writhe using exact Gauss double integral."
    global pi
    cdef int nr = r.shape[0]
    cdef double writhe_int = 0
    cdef double area, sign, dot_product
    cdef double r_df[4][3]
    cdef double r12[3]
    cdef double r34[3]
    cdef double rc[3]
    cdef double n[4][3]
    cdef np.intp_t i, j, k
    cdef bint is_zero

    for i in xrange(2, nr-1):
        for j in xrange(i-1):
            #Compute difference vectors
            for k in xrange(3):
                r_df[0][k] = r[j, k] - r[i, k]  # r13
                r_df[1][k] = r[j+1, k] - r[i, k]  # r14
                r_df[2][k] = r[j+1, k] - r[i+1, k]  # r24
                r_df[3][k] = r[j, k] - r[i+1, k]  # r23
                r12[k] = r[i+1, k] - r[i, k]
                r34[k] = r[j+1, k] - r[j, k]

            #Normal vectors
            for k in xrange(4):
                if k == 3:
                    cross(r_df[3], r_df[0], n[k])
                else:
                    cross(r_df[k], r_df[k+1], n[k])

            #Normalize the normal vectors. If any normal vectors is 0,
            #just continue.
            is_zero = False
            for k in xrange(4):
                if is_vec_nil(n[k]):
                    is_zero = True
                    break
                else:
                    normalize(n[k])
            if is_zero:
                continue

            #Compute the spherical quadrangle area
            area = 0
            for k in xrange(4):
                if k == 3:
                    dp = dot(n[3], n[0])
                else:
                    dp = dot(n[k], n[k+1])
                area += arcsin(dp)

            #Compute the sign of the area
            cross(r34, r12, rc)
            sign = dot(rc, r_df[0])
            if sign >= 0:
                writhe_int += area
            else:
                writhe_int -= area
    return writhe_int


##########################################################
# Python functions #
def writhe_exact(dr):
    '''
    writhe_exact(dr)

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
    writhe_fuller :
        Compute the writhe of the helix using the Fuller approximated single
        integral.
    ribbon_twist :
        Compute the ribbon-twist of a helix.

    Notes
    -----
    Details of the evaluation scheme is dicussed in reference [#]_ [#]_.

    References
    ----------
    .. [#] Rossetto V, Maggs AC (2003) Writhing geometry of open DNA.
        J. Chem. Phys. 118: 9864-9874.

    .. [#] Klenin K, Langowski J (2000) Computation of writhe in modeling of
        supercoiled DNA. Biopolymers 54: 307-317.
    '''
    r0 = np.vstack((np.zeros(3), np.cumsum(dr, axis=0)))  # convert to coords
    length = len(r0)
    r1 = r0[1:] - r0[0]
    r2 = r0[:(length - 1)] - r0[length - 1]
    writhe_int = _writhe_int(r0)
    return writhe_int + writhe_fuller(r1) + writhe_fuller(-r2)


def writhe_fuller(double[:, ::1] dr, bool return_val_only=True):
    '''
    writhe_fuller(dr, return_val_only=True)

    Compute the writhe using the Fuller's approximated single integral.
    Only guarantee to be correct modulo 4 pi. O(N) complexity.

    Parameters
    ----------
    dr : ndarray, shape = (N,3)
        Input delta_r vectors (dr[i] = r[i+1] - r[i]).
    return_val_only : bool, optional
        If True, return one float value for the overall writhe of the system.
        Otherwise return the individual writhe contribution for each bp-step.
        Default set to True.

    Returns
    -------
    writhe : float or ndarray of shape(N)
        The overall writhe of the helix or the individual writhe contribution
        for each bp-step. See the 'return_val_only' parameter above.

    See Also
    --------
    writhe_exact :
        Compute the writhe of the helix using the exact Gauss double integral.
    ribbon_twist :
        Compute the ribbon-twist of a helix.

    Notes
    -----
    See `writhe_exact` for references on details of the evaluation scheme.
    '''
    global pi
    cdef int n_area = dr.shape[0] - 1
    cdef double norm1, norm2, sign, dp, a
    cdef double n[3][3]
    cdef np.intp_t i, j
    cdef bint is_zero
    cdef double[::1] area = np.empty(n_area)

    for i in xrange(n_area):
        #Compute normal vectors
        cross_ez_rev(&dr[i, 0], n[0])
        cross(&dr[i, 0], &dr[i+1, 0], n[1])
        cross_ez(&dr[i+1, 0], n[2])

        #Normalize the normal vectors. If any normal vectors is 0,
        #just continue.
        is_zero = False
        for j in xrange(3):
            if is_vec_nil(n[j]):
                is_zero = True
                break
            else:
                normalize(n[j])
        if is_zero:
            area[i] = 0
            continue

        #Compute the spherical trangle area
        a = 0.5 * pi
        for j in xrange(3):
            if j == 2:
                dp = dot(n[2], n[0])
            else:
                dp = dot(n[j], n[j+1])
            a += arcsin(dp)

        #Sign of the triangle
        sign = -dr[i, 1] * dr[i+1, 0] + dr[i, 0] * dr[i+1, 1]
        if sign >= 0:
            area[i] = a
        else:
            area[i] = -a

    if return_val_only:
        return np.sum(area)
    else:
        return area


def ribbon_twist(dr, rb_vec, return_val_only=True, twist_center=0.0):
    '''
    ribbon_twist(dr, rb_vec, return_val_only=True, twist_center=0.0)

    Compute the ribbon-twist (supercoiling twist) of a helix.

    Parameters
    ----------
    dr : ndarray, shape (N,3)
        Input delta_r vectors (dr[i] = r[i+1] - r[i]).
    rb_vec : ndarray, shape (N+1,3)
        The ribbon vectors for the helix, equals to ``frame[:,:,1]``
        in the current setting. Must be normalized vectors.
    return_val_only : bool, optional
        If True, return one float value for the overall twist of the system.
        Otherwise return the individual twists for each bp-step.
        Default set to True.
    twist_center : float, optional
        Fold the step twists to (center - pi, center + pi], default to 0.

    Returns
    -------
    twist : float or ndarray of shape (N)
        The overall twist or the individual twists of the helix.
        See the 'return_val_only' parameter above.

    Raises
    ------
    ValueError
        If dr.shape[0] != rb_vec.shape[0] + 1.

    See Also
    --------
    writhe_fuller :
        Compute the writhe of the helix using the Fuller approximated single
        integral.
    writhe_exact :
        Compute the writhe of the helix using the exact Gauss double integral.
    '''
    if dr.shape[0] == rb_vec.shape[0] - 1:
        dr = np.vstack((ez, dr, ez))
    if dr.shape[0] != rb_vec.shape[0] + 1:
        raise ValueError('dr.shape[0] != rb_vec.shape[0] + 1')

    global pi
    cdef double[:, ::1] dr_c = dr
    cdef double[:, ::1] rb_vec_c = rb_vec.copy()
    cdef int n_dr = dr_c.shape[0]
    twist = np.empty(n_dr - 2)
    cdef double[::1] twist_c = twist
    cdef double b0[3]
    cdef double b1[3]
    cdef double rc[3]
    cdef double angle, sign, alpha0, alpha1, beta, twist_temp
    cdef np.intp_t i

    # Decide the range of step twist
    cdef double twist_lower = twist_center - pi
    cdef double twist_upper = twist_center + pi

    #Compute inital b vector and alpha
    cross(&dr_c[0, 0], &dr_c[1, 0], b0)
    if is_vec_nil(b0):
        normal_vec(& dr_c[0, 0], b0)
    normalize(b0)

    angle = arccos(dot(b0, &rb_vec_c[0, 0]))
    cross(b0, &rb_vec_c[0, 0], rc)
    sign = dot(rc, &dr_c[0, 0])  # sign = dot(dr[0], cross(b0, rb_vec[0]))
    if sign >= 0:
        alpha0 = angle
    else:
        alpha0 = -angle

    #Loop through the helix
    for i in xrange(n_dr-2):
        #b vector
        cross(&dr_c[i+1, 0], &dr_c[i+2, 0], b1)
        if is_vec_nil(b1):
            normal_vec(&dr_c[i+1, 0], b1)
        normalize(b1)

        #alpha angle
        angle = arccos(dot(b1, &rb_vec_c[i+1, 0]))
        cross(b1, &rb_vec_c[i+1, 0], rc)
        sign = dot(rc, &dr_c[i+1, 0])
        if sign >= 0:
            alpha1 = angle
        else:
            alpha1 = -angle

        #beta angle
        angle = arccos(dot(b0, b1))
        cross(b0, b1, rc)
        sign = dot(rc, &dr_c[i+1, 0])  # sign = dot(dr[i+1], cross(b0, b1))
        if sign >= 0:
            beta = angle
        else:
            beta = -angle

        #Compute twist
        twist_temp = beta + alpha1 - alpha0
        if twist_temp > twist_upper:
            twist_temp -= 2 * pi
        elif twist_temp <= twist_lower:
            twist_temp += 2 * pi
        twist_c[i] = twist_temp

        #Move to the next segment
        alpha0 = alpha1
        b0[0] = b1[0]
        b0[1] = b1[1]
        b0[2] = b1[2]
    if return_val_only:
        return np.sum(twist)
    else:
        return twist
##########################################################
