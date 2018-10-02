// Copyright (C) 2018 University of Vienna
// All rights reserved.
// BSD license.
// Author: Ali Baharev <ali.baharev@gmail.com>
/*
gcc -c -Ofast -march=native -std=c99 -fPIC -Wall -Wextra perturb.c
gcc -shared -Ofast perturb.o -o perturb_impl.so
rm perturb.o

------------------------------------------------------------------------
Debug build:

clang -c -ggdb3 -fsanitize=address -fno-omit-frame-pointer \
      -O0 -std=c99 -fPIC -Wall -Wextra -Wno-unused-variable \
      -Wno-unused-parameter perturb.c
clang -ggdb3 -fsanitize=address -shared -fsanitize=address -shared-libasan \
      -O0 -shared perturb.o -o perturb_impl.so
rm perturb.o

*/

static void mat_vec_mult(const double* A_flat, 
                         const int m, 
                         const int n, 
                         const double* b, 
                         double* c)
{
    const double (*A)[n] = (const double(*)[n]) A_flat;
    for (int i=0; i<m; ++i) {
        double sum = 0.0;
        for (int j=0; j<n; ++j)
            sum += A[i][j] * b[j];
        c[i] = sum;
    }
}

void perturb_impl(
        const double* __restrict__ dx_new_flat, 
        const int n_pts,
        const int m,
        const int n1,
        const int n2,
        const double* __restrict__ J33,
        const double* __restrict__ Ainv,
        const double* __restrict__ A,
        const double* __restrict__ x_subp,
        double* __restrict__ dr_norm,
        double* __restrict__ x_pert_flat
    )
{
/*  for k, dx3 in enumerate(dx_new):
        b[-n_new:] = J33 @ dx3
        dx1_dx2 = Ainv @ b
        dr = A @ dx1_dx2 - b
        dr_norm[k] = np.dot(dr, dr)
        x_pert[k] = x[x_slc.subp] + np.concatenate((-dx1_dx2, dx3)) */
    const int n = n1 + n2;
    double b[m];
    for (int i=0; i<m; ++i)
        b[i] = 0.0;
    double dx1_dx2[n1];
    double dr[m];
    const double (*dx_new)[n2] = (const double(*)[n2]) dx_new_flat;
    double (*x_pert)[n] = (double(*)[n]) x_pert_flat;
    // for k, dx3 in enumerate(dx_new):
    for (int k=0; k<n_pts; ++k) {
        const double* dx3 = dx_new[k];
        // b[-n_new:] = J33 @ dx3
        mat_vec_mult(J33, n2, n2, dx3, b + (m - n2));
        // dx1_dx2 = Ainv @ b
        mat_vec_mult(Ainv, n1, m, b, dx1_dx2);
        // dr = A @ dx1_dx2 - b
        mat_vec_mult(A, m, n1, dx1_dx2, dr);
        for (int i=m-n2; i<m; ++i)
            dr[i] -= b[i];
        // dr_norm[k] = np.dot(dr, dr)
        double sum = 0.0;
        for (int i=0; i<m; ++i)
            sum += dr[i]*dr[i];
        dr_norm[k] = sum;
        // x_pert[k] = x[x_slc.subp] + np.concatenate((-dx1_dx2, dx3))
        for (int i=0; i<n1; ++i)
            x_pert[k][i] = x_subp[i] - dx1_dx2[i];
        for (int i=n1; i<n; ++i)
            x_pert[k][i] = x_subp[i] + dx3[i-n1];
    }
}
