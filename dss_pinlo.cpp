#include <cassert>
#include <fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_interp.h>
#include "dss_pinlo.h"
#include "interp2d.h"

DSSpiNLO::DSSpiNLO(const char* filename) : filename(filename) {
    // read data from the file
    std::cerr << "Reading FF data from file " << filename << std::endl;
    double lnz;
    double lnqs2;
    size_t index;
    std::ifstream input;
    
    input.open(filename);
    
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
    }
}

void DSSpiNLO::update(double z, double qs2) {
    lnz = log(z);
    lnqs2 = log(qs2);
    // loop takes care of u, ubar, d, dbar, s, sbar, c??, cbar??, gluons
    // all fragmenting to pi+
    for (size_t i = 0; i < number_of_flavors; i++) {
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
