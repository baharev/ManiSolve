# Copyright (C) 2017 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
from __future__ import print_function, division
from cffi import FFI
from numpy import float64, intc

__all__ = ['CProblem']

ffi = FFI()

ffi.cdef('''

void evaluate(int index, const double* py_point, double* py_residual);

void jacobian_evaluation(int index, 
                         const double* py_point, 
                         double* py_residual, 
                         double* py_2d_array);

void solve(int index,
           double* x_2d_array, 
           double* res_2d_array,
           const int n_points,
           const int* idx_fixed_2d,
           const double* val_fixed_2d,
           const int n_cols_fixed_2d,
           const double tolerance, // 1.0e-8, >= EPS[i]^2
           const int n_trials,
           unsigned int seed,
           const int iprint);

''')


NULL = ffi.NULL

def double_ptr(arr):
    assert arr.dtype == float64, arr.dtype
    return ffi.cast('double*', arr.ctypes.data)

def const_int_ptr(arr):
    assert arr.dtype == intc, arr.dtype
    return ffi.cast('const int*', arr.ctypes.data)

def const_double_ptr(arr):
    assert arr.dtype == float64, arr.dtype
    return ffi.cast('const double*', arr.ctypes.data)

def assert_2d_arrays_of_equal_length(a, b):
    ma, na = a.shape  # a must be a 2D array 
    mb, nb = b.shape  # b must be a 2D array
    assert ma == mb and na and nb, (a.shape, b.shape)


class CProblem:
    
    __slots__ = ['so']
    
    def __init__(self, so_name):
        self.so = ffi.dlopen(so_name)
    
    def evaluate(self, index, x, r):
        assert len(x.shape) == 1 and len(r.shape) == 1, '1D arrays expected'
        x = const_double_ptr(x)
        r = double_ptr(r)
        self.so.evaluate(index, x, r)
    
    def jacobian_evaluation(self, index, x, r, jac_2d_arr):
        # jac_2d_array is a NumPy array storing the Jacobian as a DENSE matrix
        m, n = jac_2d_arr.shape
        assert m <= len(r) and n <= len(x)
        x = const_double_ptr(x)
        r = double_ptr(r)
        jac_2d_arr = double_ptr(jac_2d_arr)
        self.so.jacobian_evaluation(index, x, r, jac_2d_arr)
    
    def solve(self, index, x, r, tol=1.0e-8, n_trials=10, seed=1, iprint=0):
        # The x and r are NOT touched if the solver fails, otherwise their 
        # corresponding entries will be set to the solution.
        # The x and r are 2D arrays: set of points and the corresponding residuals.
        # The seed is not set (changed) if seed is equal to 0.
        # The tolerance should be at least EPS[i]^2 (=1.0e-8), see the FORTRAN code.
        assert_2d_arrays_of_equal_length(x, r)
        n_points = x.shape[0]
        x = double_ptr(x)
        r = double_ptr(r)
        self.so.solve(index, x, r, n_points, NULL, NULL, 0, tol, n_trials, seed, iprint)
    
    def solve_from(self, index, x, r, tol=1.0e-8, n_trials=1, seed=1, iprint=0):
        self.solve(index, x, r, tol, -n_trials, seed, iprint) # note the flipped sign
    
    def solve_fixed(self, index, x, r, idx, val, tol=1.0e-8, n_trials=40, seed=1, iprint=0):
        # See the comments just above at solve()!
        # Fixed variables: for point i, idx[i] gives their index array (index in the
        # va27_point), val[i] their value array.
        assert_2d_arrays_of_equal_length(x, r)
        assert_2d_arrays_of_equal_length(idx, val)
        assert x.shape[0] == idx.shape[0], (x.shape, idx.shape)    
        n_points = x.shape[0]
        n_fixed = val.shape[1]
        x = double_ptr(x)
        r = double_ptr(r)
        idx = const_int_ptr(idx)
        val = const_double_ptr(val)
        self.so.solve(index, x, r, n_points, idx, val, n_fixed, tol, n_trials, seed, iprint)
    
    def solve_fixed_from(self, index, x, r, idx, val, tol=1.0e-8, n_trials=1, seed=1, iprint=0):
        self.solve_fixed(index, x, r, idx, val, tol, -n_trials, seed, iprint) # note the flipped sign
