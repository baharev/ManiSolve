// Copyright (C) 2017 University of Vienna
// All rights reserved.
// BSD license.
// Author: Ali Baharev <ali.baharev@gmail.com>
#include <assert.h>
#include <float.h>
#include <math.h>

/* Compile with:


Consider adding -DNDEBUG for further speed:

gcc -c -Ofast -march=native -std=c99 -fPIC -Wall -Wunused-variable -Wextra subsampling2.c
gcc -shared -Ofast subsampling2.o -o subsample2.so
rm subsampling2.o

-----------------------------------------------------------------------

This does not make much sense since there is no allocation in this code:

gcc -c -ggdb3 -fsanitize=address -fno-omit-frame-pointer -O0 -std=c99 -fPIC -Wall -Wextra subsampling.c
gcc -shared -O0 subsampling.o -o subsample.so -Wl,--whole-archive /usr/lib/gcc/x86_64-linux-gnu/5/libasan.a -Wl,--no-whole-archive
rm subsampling.o

No longer required with gcc 4.9.3 and above:
python subsampling.py 2>&1 | misc/va27/asan_symbolize.py

*/

static int max_dmin(const double* d_min,
                    const int n,
                    const unsigned char* selected)
{
    // Find the max element in the array of the nearest neighbors for the
    // not yet selected points.
    int i_max = -1;
    double d_max = -1.0;

    for (int i=0; i<n; ++i) {
        if (!selected[i]) {
            double d = d_min[i];
            if (d > d_max) {
                i_max = i;
                d_max = d;
            }
        }
    }
    assert(i_max != -1);
    assert(d_min[i_max] >= 0.0);

    return i_max;
}
/*
static double distance(const double* x, int point_dim, int i, int j) {
    // Distance in the infinity norm
    const double* a = x + i*point_dim;
    const double* b = x + j*point_dim;
    double dist = 0.0;
    for (int k=0; k<point_dim; ++k) {
        const double tmp = fabs(a[k] - b[k]);
        if (tmp > dist)
            dist = tmp;
    }
    return dist;
}
*/
static double distance(const double* x, int point_dim, int i, int j) {
    // Distance in the L1-norm
    const double* a = x + i*point_dim;
    const double* b = x + j*point_dim;
    double dist = 0.0;
    for (int k=0; k<point_dim; ++k) {
        dist += fabs(a[k] - b[k]);
    }
    return dist;
}

static void add_points(const double* __restrict x,
                       unsigned char* __restrict__ selected,
                       const int n,
                       const int point_dim,
                       double* __restrict__ d_min)
{
    // Add that non-selected point next whose nearest neighbor among the
    // selected points is the most distant.
    const int index = max_dmin(d_min, n, selected);

    assert(index != -1);
    assert(!selected[index]);

    selected[index] = 1;

    // Update the nearest neighbors. Is the newly selected one the nearest now?
    for (int i=0; i<n; ++i) {
        if (!selected[i]) { // <-- Can become a bottleneck. Track the selected indices instead then?
            const double d = distance(x, point_dim, i, index);
            if (d < d_min[i]) {
                d_min[i] = d;
            }
        }
    }
}

static void set_d_to_nn(const double* __restrict x,               // shape: n x point_dim
                        const int* __restrict__ idx_selected,     // shape: k
                        const int k,
                        const int* __restrict__ idx_not_selected, // shape: n-k
                        const int n,                              // number of points
                        const int point_dim,
                        double* __restrict__ d_min)
{
    // Initialize the array of nearest neighbors
    for (int ii=0; ii < (n-k); ++ii) {
        // for each not selected point
        const int i = idx_not_selected[ii];
        // find its nearest neighbor among the selected ones
        double val_min = DBL_MAX;
        for (int jj=0; jj < k; ++jj) {
            const int j = idx_selected[jj];
            double d = distance(x, point_dim, i, j);
            if (d < val_min)
                val_min = d;
        }
        assert(val_min >= 0.0);
        assert(val_min < DBL_MAX);
        d_min[i] = val_min;
    }
}

void downsample2(const double* __restrict__ x,             // shape: n x point_dim
                 const int* __restrict__ idx_selected,     // shape: k
                 int k,
                 const int* __restrict__ idx_not_selected, // shape: n-k
                 unsigned char* __restrict__ selected,     // shape: n
                 double* __restrict__ d_min,               // shape: n
                 int n,                                    // number of points
                 int point_dim,
                 int n_pts_to_add)
{
    // If the calling Python code accidentally passes in a pointer to a view,
    // we are reading the wrong data. Try to detect it on a best effort basis.
    for (int i=0; i<n*point_dim; ++i)
        assert(isfinite(x[i]) && 
        "NaN or inf among the point coordinates; x must be contiguous, and not a view.");
    
    // See the calling Python code for the preconditions, they are checked there
    set_d_to_nn(x, idx_selected, k, idx_not_selected, n, point_dim, d_min);

    for (int i=0; i<n_pts_to_add; ++i) {
        add_points(x, selected, n, point_dim, d_min);
    }
}
