# Copyright (C) 2018 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from cffi import FFI
from numpy import float64, ascontiguousarray

ffi = FFI()

ffi.cdef('''

void perturb_impl(
        const double* dx_new_flat, 
        const int n_pts,
        const int m,
        const int n1,
        const int n2,
        const double* J33,
        const double* Ainv,
        const double* A,
        const double* x_subp,
        double* dr_norm,
        double* x_pert_flat
    );

''')


NULL = ffi.NULL

so = ffi.dlopen('perturb_impl.so')
perturb_impl = so.perturb_impl

def double_ptr(arr):
    assert arr.dtype == float64, arr.dtype
    return ffi.cast('double*', arr.ctypes.data)

def const_double_ptr(arr):
    assert arr.dtype == float64, arr.dtype
    return ffi.cast('const double*', arr.ctypes.data)

#-------------------------------------------------------------------------------
#     for k, dx3 in enumerate(dx_new):
#         b[-n_new:] = J33 @ dx3
#         dx1_dx2 = Ainv @ b
#         dr = A @ dx1_dx2 - b
#         dr_norm[k] = np.dot(dr, dr)
#         x_pert[k] = x[x_slc.subp] + np.concatenate((-dx1_dx2, dx3))

def perturb_C(dx_new, J33, Ainv, A, x_subp, dr_norm, x_pert):
    J33 = ascontiguousarray(J33)
    A = ascontiguousarray(A)
    assert dx_new.flags['C_CONTIGUOUS']  == True, dx_new.flags
    assert J33.flags['C_CONTIGUOUS']     == True, J33.flags
    assert Ainv.flags['C_CONTIGUOUS']    == True, Ainv.flags
    assert A.flags['C_CONTIGUOUS']       == True, A.flags
    assert x_subp.flags['C_CONTIGUOUS']  == True, x_subp.flags
    assert dr_norm.flags['C_CONTIGUOUS'] == True, dr_norm.flags
    assert x_pert.flags['C_CONTIGUOUS']  == True, x_pert.flags
    m, n1 = A.shape 
    assert Ainv.shape == (n1, m), (Ainv.shape, m, n1)
    assert J33.shape[0] == J33.shape[1], J33.shape
    n2 = J33.shape[0]
    assert x_subp.ndim == 1, x_subp.shape
    n = len(x_subp)
    assert n == n1 + n2, (n, n1, n2)
    n_pts, n_new = dx_new.shape
    assert n_pts == len(dr_norm), (n_pts, len(dr_norm))
    assert n_pts == len(x_pert), (n_pts, len(x_pert))
    assert n_new == n2, (n_new, n2)
    assert x_pert.shape[1] == n
    perturb_impl(
        const_double_ptr(dx_new), 
        n_pts,
        m,
        n1,
        n2,
        const_double_ptr(J33),
        const_double_ptr(Ainv),
        const_double_ptr(A),
        double_ptr(x_subp),
        double_ptr(dr_norm),
        double_ptr(x_pert)
    )
