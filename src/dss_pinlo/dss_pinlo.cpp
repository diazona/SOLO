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

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <gsl/gsl_interp.h>
#include "dss_pinlo.h"
#include "interp2d.h"

DSSpiNLO::DSSpiNLO(const char* filename) :
 number_of_lnz_values(0), number_of_lnqs2_values(0),
 lnz_array(NULL), lnqs2_array(NULL), ff_arrays(NULL),
 filename(filename) {
    std::cerr << "Reading FF data from file " << filename << std::endl;
    double lnz;
    double lnqs2;
    size_t index;
    std::ifstream input;
    
    size_t lines = 0;
    std::set<double> lnz_values;
    std::set<double> lnqs2_values;
    
    input.open(filename);
    if (!input.good()) {
        throw std::ios_base::failure("Unable to read file");
    }
    
    // Read through once to identify unique lnz and lnQs2 values
    input >> lnz >> lnqs2;
    while (!input.eof()) {
        input.ignore(256, '\n'); // read and ignore until end of line
        lines++;
        lnz_values.insert(lnz);
        lnqs2_values.insert(lnqs2);
        input >> lnz >> lnqs2;
    }
    number_of_lnz_values = lnz_values.size();
    number_of_lnqs2_values = lnqs2_values.size();
    if (lines != number_of_lnz_values * number_of_lnqs2_values) {
        input.close();
        std::cerr << "Data are not arranged in a complete grid" << std::endl;
        std::cerr << "Found " << number_of_lnz_values << " values for ln(z)," << std::endl;
        std::cerr << "      " << number_of_lnqs2_values << " values for ln(Qs2)," << std::endl;
        std::cerr << "      " << lines << " lines read from file" << std::endl;
        exit(1);
    }
    
    lnz_array = new double[number_of_lnz_values];
    lnqs2_array = new double[number_of_lnqs2_values];
    ff_arrays = new double*[number_of_flavors];
    for (size_t i = 0; i < number_of_flavors; i++) {
        ff_arrays[i] = new double[number_of_lnz_values * number_of_lnqs2_values];
    }
    
    input.clear(); // have to clear the EOF bit before seeking
    input.seekg(0, std::ios::beg);
    
    // Read through again to actually get the data
    for (size_t iz = 0; iz < number_of_lnz_values; iz++) {
        for (size_t iq = 0; iq < number_of_lnqs2_values; iq++) {
            index = INDEX_2D(iz, iq, number_of_lnz_values, number_of_lnqs2_values);
            input >> lnz >> lnqs2;
            input >> ff_arrays[up][index]
                  >> ff_arrays[up_bar][index]
                  >> ff_arrays[down][index]
                  >> ff_arrays[down_bar][index]
                  >> ff_arrays[strange][index]
                  >> ff_arrays[strange_bar][index]
                  >> ff_arrays[charm][index]
                  >> ff_arrays[charm_bar][index]
                  >> ff_arrays[gluon][index];
            if (iz == 0) {
                lnqs2_array[iq] = lnqs2;
            }
        }
        lnz_array[iz] = lnz;
    }
    input.close();
    
    std::cerr << "Constructing interpolation" << std::endl;
    
    // construct interpolation objects
    for (size_t i = 0; i < number_of_flavors; i++) {
        lnz_accel[i] = gsl_interp_accel_alloc();
        lnqs2_accel[i] = gsl_interp_accel_alloc();
        // TODO: change to bicubic once it's tested
        interpolators[i] = interp2d_alloc(interp2d_bilinear, number_of_lnz_values, number_of_lnqs2_values);
        interp2d_init(interpolators[i], lnz_array, lnqs2_array, ff_arrays[i], number_of_lnz_values, number_of_lnqs2_values);
    }
    
    std::cerr << "Done initializing DSSpiNLO" << std::endl;
}

DSSpiNLO::~DSSpiNLO() {
    for (size_t i = 0; i < number_of_flavors; i++) {
        gsl_interp_accel_free(lnz_accel[i]);
        lnz_accel[i] = NULL;
        gsl_interp_accel_free(lnqs2_accel[i]);
        lnqs2_accel[i] = NULL;
        interp2d_free(interpolators[i]);
        interpolators[i] = NULL;
        delete[] ff_arrays[i];
        ff_arrays[i] = NULL;
    }
    delete[] ff_arrays;
    ff_arrays = NULL;
    delete[] lnz_array;
    lnz_array = NULL;
    delete[] lnqs2_array;
    lnqs2_array = NULL;
}

void DSSpiNLO::update(double z, double qs2) {
    lnz = log(z);
    lnqs2 = log(qs2);
    // loop takes care of u, ubar, d, dbar, s, sbar, c??, cbar??, gluons
    // all fragmenting to pi+
    for (size_t i = 0; i < number_of_flavors; i++) {
        if (lnz < lnz_array[0] || lnz > lnz_array[number_of_lnz_values-1] || lnqs2 < lnqs2_array[0] || lnqs2 > lnqs2_array[number_of_lnqs2_values-1]) {
            throw FragmentationFunctionRangeException(z, qs2);
        }
        pi_plus_ff[i] = interp2d_eval(interpolators[i], lnz_array, lnqs2_array, ff_arrays[i], lnz, lnqs2, lnz_accel[i], lnqs2_accel[i]) / z;
    }
    pi_minus_ff[up] = pi_plus_ff[up_bar];
    pi_minus_ff[up_bar] = pi_plus_ff[up];
    pi_minus_ff[down] = pi_plus_ff[down_bar];
    pi_minus_ff[down_bar] = pi_plus_ff[down];
    pi_minus_ff[strange] = pi_plus_ff[strange_bar];
    pi_minus_ff[strange_bar] = pi_plus_ff[strange];
    pi_minus_ff[charm] = pi_plus_ff[charm_bar];
    pi_minus_ff[charm_bar] = pi_plus_ff[charm];
    pi_minus_ff[gluon] = pi_plus_ff[gluon];
    for (size_t i = 0; i < number_of_flavors; i++) {
        pi_zero_ff[i] = 0.5 * (pi_plus_ff[i] + pi_minus_ff[i]);
    }
}

// charge is +1, -1, 0 to indicate pion charge
double DSSpiNLO::fragmentation(flavor f, hadron h) {
    switch(h) {
        case pi_plus:
            return pi_plus_ff[f];
        case pi_zero:
            return pi_zero_ff[f];
        case pi_minus:
            return pi_minus_ff[f];
        default:
            assert(false);
    }
}
