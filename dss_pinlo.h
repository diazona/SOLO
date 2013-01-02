/*
 * A C++ interface to the DSS pion fragmentation functions.
 * 
 * Copyright 2012 David Zaslavsky
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

#ifndef DSS_PINLO_H_INCLUDED
#define DSS_PINLO_H_INCLUDED

#include "interp2d.h"

class DSSpiNLO {
private:
    const static size_t number_of_flavors = 9;

    size_t number_of_lnz_values;
    size_t number_of_lnqs2_values;

    double* lnz_array;
    double* lnqs2_array;
    double** ff_arrays;

    const char* filename;
    
    interp2d* interpolators[number_of_flavors];
    gsl_interp_accel* lnz_accel[number_of_flavors];
    gsl_interp_accel* lnqs2_accel[number_of_flavors];

    double lnz;
    double lnqs2;
    double pi_plus_ff[number_of_flavors];
    double pi_minus_ff[number_of_flavors];
    double pi_zero_ff[number_of_flavors];
    
public:
    enum flavor {gluon, up, up_bar, down, down_bar, strange, strange_bar, charm, charm_bar};
    enum hadron {pi_plus, pi_zero, pi_minus};

    DSSpiNLO(const char* filename);
    ~DSSpiNLO();
    void update(double z, double qs2);
    double fragmentation(flavor f, hadron h);
};

#endif