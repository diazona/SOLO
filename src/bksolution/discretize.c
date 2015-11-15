/*
 * Part of SOLO
 *
 * Copyright 2015 David Zaslavsky
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <math.h>
#include <gsl/gsl_integration.h>

/**
 * A quick check whether the linear basis function is equal to zero.
 *
 * This exists because it's useful for speed to check whether the
 * function will be equal to zero first, so we can avoid evaluating
 * either the linear basis function, or the function being integrated.
 *
 * @param[in] x the value at which to evaluate T_n
 * @param[in] n the index n (of course)
 * @param[in] step the step δx
 */
inline int linear_basis_function_is_zero(const double x, const int n, const double step) {
    return (x < (p->n - 1) * p->step) || (x > (p->n + 1) * p->step);
}

/**
 * A basis function for linear interpolation.
 *
 * Suppose there is a function f which is sampled at regularly spaced grid
 * points nδx, and we have f_{n} = f(nδx) for each n. We want to be able to
 * approximate f(x) for any x via linear interpolation. The typical way to
 * do this, given a value of x, is to find the index m such that
 * mδx < x < (m+1)δx, and calculate the linear combination
 *  f(x) ~ f_{m} + (x - mδx)(f_{m+1} - f_{m})/(δx)
 * But it's also possible (and mathematically equivalent) to use basis
 * functions T_n(x) to express
 *  f(x) ~ sum_n f_{n} T_n(x)
 * For linear interpolation, these basis functions are triangular, peaked
 * at x_n and equal to zero outside [x_{n-1},x_{n+1}], so this works out to
 *  f(x) ~ f_{m} T_m(x) + f_{m+1} T_{m+1}(x)
 * A little algebra shows that
 *  T_m(x) = 1 - (mδx - x)/δx,  (m-1)δx < x < mδx
 *  T_m(x) = 1 - (x - mδx)/δx,  mδx < x < (m+1)δx
 *  T_m(x) = 0,                 otherwise
 * This function implements T_n(x).
 *
 * @param[in] x the value at which to evaluate T_n
 * @param[in] n the index n (of course)
 * @param[in] step the step δx
 */
inline double linear_basis_function(const double x, const int n, const double step) {
    if (linear_basis_function_is_zero(x, n, step)) {
        return 0;
    }
    double peak = n * step;
    return 1 - fabs(x - peak)/step;
}

struct parameters {
    int n;
    double step;

    double (*integrand)(double, void*);
    void* integrand_params;
};

/**
 * Implements the integrand for linear discretization, a product of
 * the function being integrated and a linear basis function.
 */
double linear_discretization_integrand(const double x, void* params) {
    struct parameters* p = (struct parameters*)params;

    if (linear_basis_function_is_zero(x, p->n, p->step)) {
        return 0;
    }
    else {
        return linear_basis_function(x, p->n, p->step) * p->integrand(x, p->integrand_params);
    }
}

void linear_discretize(double (*function)(double, void*), void* params, const double start, const double step, const size_t n_points, const size_t subinterval_limit, double* result, double* abserr) {
    if (result == NULL) {
        return;
    }
    double l_res, l_err;
    struct parameters discretization_params = {0, step, function, params};
    gsl_function func;
    func.function = linear_discretization_integrand;
    func.params = &discretization_params;
    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(subinterval_limit);
    double istep = 1.0 / step;

    for (size_t n = 0; n < n_points; n++) {
        int status;
        discretization_params.n = n;
        status = gsl_integration_qag(func, start + (n - 1) * step, start + (n + 1) * step, 0, 1e-5, subinterval_limit, GSL_INTEG_GAUSS61, workspace, &l_res, &l_err);
        result[n] = (status == GSL_SUCCESS) ? l_res * istep : NAN;
        if (abserr != NULL) {
            abserr[n] = (status == GSL_SUCCESS) ? l_err * istep : NAN;
        }
    }

    gsl_integration_workspace_free(workspace);
}
