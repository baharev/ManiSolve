#include <assert.h>
#include <math.h>
#include <stdlib.h>

static double x_guess(double lb, double ub, double r) {
    assert(lb < ub);
    assert(lb >= 0.0);
    assert(0.0 <= r && r <= 1.0);
    const double m = ub - lb + 1.0;
    return pow(m, r) + lb - 1.0;
}

static double x_initial(double lb, double ub, double r) {

    if (lb >= 0.0)
        return  x_guess(lb, ub, r);
    if (ub <= 0.0)
        return -x_guess(-ub, -lb, r);
    assert(lb < 0.0 && 0.0 < ub);
    double a = log(-lb + 1.0);
    double b = log( ub + 1.0);
    double s = a / (a + b);
    if (r <= s)
        return -x_guess(0.0, -lb, (s-r)/s);
    else
        return x_guess(0.0, ub, (r-s)/(1.0-s));    
}


void randpoint(const double* __restrict__ lb, 
               const double* __restrict__ ub, 
               const int n, 
                     double* __restrict__ x) 
{
    const double scale = 1.0 / (double) RAND_MAX;
    
    for (int i=0; i<n; ++i)
        x[i] = x_initial(lb[i], ub[i], rand()*scale);
}

void set_seed(unsigned int seed) {
    srand(seed);   
}
