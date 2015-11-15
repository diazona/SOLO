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

#include <cmath>
#include <gsl/gsl_math.h>
#include "solver.h"

#define ln(x) log(x)

using std::log;

double splitting(double Y, int alpha, void* params) {
}

double kernel(double k2, double kprime2, int alpha, void* params) {
    switch (alpha) {
        case 0:
            return (kprime2 - k2) / fabs(kprime2 - k2) + k2 / sqrt(4 * gsl_pow_2(kprime2) + gsl_pow_2(k2));
        default:
            return 0;
    }
}

double delta_initial_condition(double k2, double Y, void* params) {
    if (k2 == 1.0) {
        double step = *(static_cast<double*>(params));

    }
}

int main(const int argc, const char** argv) {
    Solver s(
        delta_initial_condition, NULL,
        kernel, NULL,
        splitting, NULL,
        ln(1e-3), ln(1e3), ln(10)/10,
        -ln(1e-2), -ln(1e-5), ln(10)/10,
        1);
}
