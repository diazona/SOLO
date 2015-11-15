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

#include <stdlib.h>

/**
 * Convolve the given function with the linear basis functions (T_n) for each
 * value of n from 0. After returning, result[n] contains the result of
 * integrating T_n(x)*function(x).
 *
 * This implements the discretization procedure for linear interpolation.
 * For example, equation (73) in hep-ph/0307188.
 *
 * @param[in] function the function to be discretized
 * @param[in] params an arbitrary object to be passed to the function
 * @param[in] start the lowest value of the indpendent variable at which
 * to produce a discretized result
 * @param[in] step the step between values of the independent variable at
 * which to produce discretized results
 * @param[in] n_points the number of values at which to produce discretized
 * results
 * @param[in] subinterval_limit the maximum number of subintervals to be used
 * by the adaptive integration
 * @param[out] result an array in which to store the results, which must have
 * length at least `n_points`
 * @param[out] abserr an array in which to store the error estimates, which
 * must have length at least `n_points`, though this can be `NULL` in which
 * case the error estimates will be discarded
 */
void linear_discretize(double (*function)(double, void*), void* params, const double start, const double step, const size_t n_points, const size_t subinterval_limit, double* result, double* abserr);
