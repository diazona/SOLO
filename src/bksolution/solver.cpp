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

/* This implements the iterative relaxation method described in
 * arxiv:hep-ph/9702418, section 3, and arxiv:hep-ph/0307188, section 4.2.
 */

#include <algorithm>
#include <gsl/gsl_integration.h>
#include "../configuration/context.h"
#include "discretize.h"
#include "solver.h"

using std::copy;
using std::begin;
using std::end;

Solver::Solver(
    double (*initial_condition)(double, double, void*), void* initial_condition_params,
    double (*kernel)(double, double, int, void*), void* kernel_params,
    double (*splitting)(double, int, void*), void* splitting_params,
    double lnk_min, double lnk_max, double lnk_step,
    double Y_min, double Y_max, double Y_step,
    size_t kernel_terms,
    double subinterval_limit = 10000
) :
    initial_condition(initial_condition), initial_condition_params(initial_condition_params),
    kernel(kernel), kernel_params(kernel_params),
    splitting(splitting), splitting_params(splitting_params),
    lnk_min(lnk_min), lnk_max(lnk_max), lnk_step(lnk_step), lnk_dimension(0),
    Y_min(Y_min), Y_max(Y_max), Y_step(Y_step), Y_dimension(0),
    kernel_terms(kernel_terms),
    subinterval_limit(subinterval_limit)
{
    lnk_dimension = ceil((lnk_max - lnk_min) / lnk_step);
    assert(lnk_min + lnk_dimension * lnk_step >= lnk_max);
    Y_dimension = ceil((Y_max - Y_min) / Y_step);
    assert(Y_min + Y_dimension * Y_step >= Y_max);
}

// avoid dealing with the cruft of C++ classes - not sure it makes a difference
extern "C" {

struct initial_condition_params {
    double (*initial_condition)(double, double, void*);
    double Y;
    void* params;
};

double initial_condition_wrapper(double k2, void* params) {
    struct initial_condition_params* icp = static_cast<initial_condition_params*>(params);
    return icp->initial_condition(k2, icp->Y, icp->params);
}

struct kernel_params {
    double (*kernel)(double, double, int, void*);
    double k2;
    int alpha;
    void* params;
};

double kernel_wrapper(double kprime2, void* params) {
    struct kernel_params* kp = static_cast<kernel_params*>(params);
    return kp->kernel(kp->k2, kprime2, kp->alpha, kp->params);
}

struct splitting_params {
    double (*splitting)(double, int, void*);
    int alpha;
    void* params;
};

double splitting_wrapper(double Y, int alpha, void* params) {
    struct splitting_params* sp = static_cast<splitting_params*>(params);
    return sp->splitting(Y, sp->alpha, sp->params);
}

} // extern "C"

void Solver::setup() {
    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(subinterval_limit);

    // kernel K
    kernel_array = new double[kernel_terms][lnk_dimension][lnk_dimension];
    {
        struct kernel_params params = {kernel, 0, 0, kernel_params};
        for (size_t alpha = 0; alpha < kernel_terms; alpha++) {
            params.alpha = alpha;
            for (size_t m = 0; m < lnk_dimension; m++) {
                params.k2 = exp(2 * (lnk_min + m * lnk_step));
                linear_discretize(kernel_wrapper, params, 2 * lnk_min, 2 * lnk_step, lnk_dimension, subinterval_limit, kernel_array[m], NULL);
            }
        }
    }

    // splitting function P
    splitting_array = new double[kernel_terms][Y_dimension];
    {
        struct splitting_params params = {splitting, 0, splitting_params};
        for (size_t alpha = 0; alpha < kernel_terms; alpha++) {
            params.alpha = alpha;
            linear_discretize(splitting, splitting_params, Y_min, Y_step, Y_dimension, subinterval_limit, splitting_array, NULL);
        }
    }

    // gluon Green's function, a.k.a. BK solution, G
    solution_array[][] = new double[Y_dimension][lnk_dimension];

    // set the initial condition
    linear_discretize(initial_condition, initial_condition_params, lnk_min, lnk_step, lnk_dimension, subinterval_limit, solution_array[0], NULL);


    // run the evolution - no integrals here, this is just matrix multiplication
    for (size_t n = 0; n < Y_dimension; n++) {
        /* At each value of n, the first time through, solution_array[m][n]
         * is uninitialized. But the calculation of solution_array[m][n] (below)
         * will actually use the values of solution_array[m][n]. So we make
         * a first approximation by copying the results from the previous
         * iteration.
         */
        copy(begin(solution_array[n - 1]), end(solution_array[n - 1]), begin(solution_array[n]));

        /* Now we can use the integral equation, eq 70, to refine the
         * approximation. We keep track of how much the result changes from
         * one iteration to the next, and stop when the change becomes
         * sufficiently small.
         */
        double total_change = 0;

        do {
            for (size_t m = 0; m < lnk_dimension; m++) {
                // value of the double integral on the right side of eq 70, as approximated by eq 74
                double result = solution_array[0][m];
                for (size_t i = 0; i < lnk_dimension; i++) {
                    for (size_t j = 0; j < n - abs(m - i); j++) {
                        // I think there is no way to pass this off to the BLAS,
                        // since it's not really a matrix multiplication
                        result += splitting_array[n - j - abs(m - i)] * kernel_array[m][i] * solution_array[j][i];
                    }
                }
                /* The formula below defines the total change as follows:
                 * for each value of k (indexed by m), take the difference between
                 * the old and new results, divide it by the new result, and take the
                 * absolute value. Then add these results for all values of k. I do
                 * it this way just for simplicity. Perhaps a more sophisticated
                 * definition, like a sum of squares, would be more appropriate.
                 */
                total_change += abs((solution_array[n][m] - result) / result);

                solution_array[n][m] = result;
            }
        // why 1e-3? no reason - this is another thing that should be more carefully thought out
        } while (total_change < 1e-3 * lnk_dimension);
    }
}
