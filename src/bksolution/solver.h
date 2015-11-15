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

#include <cstdlib>

class Solver {
public:
    Solver(
        double (*initial_condition)(double, void*), void* initial_condition_params,
        double (*kernel)(double, double, void*), void* kernel_params,
        double (*splitting)(double, void*), void* splitting_params,
        double lnk_min, double lnk_max, double lnk_step,
        double Y_min, double Y_max, double Y_step,
        size_t kernel_terms,
        double subinterval_limit = 10000);

    double solution_array[][];
private:
    void setup();

    double subinterval_limit;
    double lnk_min, lnk_max, lnk_step;
    size_t lnk_dimension;
    double Y_min, Y_max, Y_step;
    size_t Y_dimension;

    double kernel_array[][][];
    double splitting_array[][];

    double (*initial_condition)(double, void*);
    void* initial_condition_params;
    double (*kernel)(double, double, void*);
    void* kernel_params;
    double (*splitting)(double, void*);
    void* splitting_params;

    size_t kernel_terms;
};
