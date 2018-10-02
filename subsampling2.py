# Copyright (C) 2017, 2018 University of Vienna
# All rights reserved.
# BSD license.
# Author: Ali Baharev <ali.baharev@gmail.com>
import numpy as np
from cffi import FFI

__all__ = [ 'subsample2' ]

ffi = FFI()
ffi.cdef('''

void downsample2(const double* x,             // shape: n x point_dim
                 const int* idx_selected,     // shape: k
                 int k,
                 const int* idx_not_selected, // shape: n-k
                 unsigned char* selected,     // shape: n
                 double* d_min,               // shape: n
                 int n,                       // number of points
                 int point_dim,
                 int n_pts_to_add);

''')

so = ffi.dlopen('subsample2.so')
downsample = so.downsample2

def const_double_ptr(arr):
    assert arr.dtype == np.float64, arr.dtype
    return ffi.cast('const double*', arr.ctypes.data)

def double_ptr(arr):
    assert arr.dtype == np.float64, arr.dtype
    return ffi.cast('double*', arr.ctypes.data)

def uchar_ptr(arr):
    assert arr.dtype == np.bool_, arr.dtype
    return ffi.cast('unsigned char*', arr.ctypes.data)

def const_int_ptr(arr):
    assert arr.dtype == np.intc, arr.dtype
    return ffi.cast('const int*', arr.ctypes.data)

#-------------------------------------------------------------------------------

def subsample2(x, n_desired_pts):
    # Returns: a boolean array, indicating whether a point in x was selected.
    assert x.ndim == 2, x.shape # x is a 2D array of the points
    n_pts, point_dim = x.shape[0], x.shape[1]    
    assert n_pts > 0 and point_dim > 0, (n_pts, point_dim)
    assert np.isfinite(x).all()
    assert 0 < n_desired_pts, n_desired_pts
    if n_desired_pts > n_pts:
        msg = '*** Warning: We have less points than desired! ({} < {}) ***'
        print(msg.format(n_pts, n_desired_pts))
        return np.full(n_pts, np.True_, dtype=np.bool_)
    else:
        # Best effort: try to pick the point, closest to the middle:
        mean = np.mean(x, axis=0)
        idx = np.argmin(np.linalg.norm(x-mean, axis=1, ord=1))
        selected = np.full(n_pts, np.False_, dtype=np.bool_)
        selected[idx] = np.True_
        # Then call subsampling to fill up the rest, if any
        return subsample2_(x, selected, n_desired_pts)

def subsample2_(x, selected, n_desired_pts):
    # Similar to subsample but some points are already selected 
    assert x.ndim == 2, x.shape # x is a 2D array of the points
    n_pts, point_dim = x.shape[0], x.shape[1]    
    assert n_pts > 0 and point_dim > 0, (n_pts, point_dim)
    assert np.isfinite(x).all()
    assert x.shape[0] == selected.shape[0], (x.shape, selected.shape)
    assert 0 < n_desired_pts, n_desired_pts
    if n_desired_pts > n_pts:
        msg = '*** Warning: We have less points than desired! ({} < {}) ***'
        print(msg.format(n_pts, n_desired_pts))
        return np.full(n_pts, np.True_, dtype=np.bool_)
    selected = selected.copy()
    all_idx = np.arange(n_pts)
    idx = all_idx[selected]
    k = len(idx)
    to_add = n_desired_pts - k
    if to_add <= 0:
        return selected
    idx_not_selected = all_idx[~selected]
    assert len(idx_not_selected) == n_pts - k
    assert selected.any()
    x = np.ascontiguousarray(x)
    idx = idx.astype(np.intc)
    idx_not_selected = idx_not_selected.astype(np.intc)
    # d_nn: temporary work array for distances to the nearest neighbors
    d_nn = np.full(n_pts, np.nan)
    downsample(const_double_ptr(x), # shape: (n_pts, point_dim)
               const_int_ptr(idx),  # shape: (k,)
               k,
               const_int_ptr(idx_not_selected), # shape: (n_pts-k,)
               uchar_ptr(selected), # shape: (n_pts,)
               double_ptr(d_nn),    # shape: (n_pts,)
               n_pts, 
               point_dim,
               to_add)
    return selected

#-------------------------------------------------------------------------------
# For performance testing and profiling with perf

def _main():
    np.random.seed(1)
    # 10k -> 1k downsampling, point dimension: 25
    subsample2(np.random.random((10000, 25)), 1000)

if __name__ == '__main__':
    _main()
