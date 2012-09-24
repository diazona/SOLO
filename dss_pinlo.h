#ifndef DSS_PINLO_H_INCLUDED
#define DSS_PINLO_H_INCLUDED

#include "interp2d.h"

class DSSpiNLO {
private:
    const static size_t number_of_lnz_values = 35;
    const static size_t number_of_lnqs2_values = 24;
    const static size_t number_of_flavors = 9;

    double lnz_array[number_of_lnz_values];
    double lnqs2_array[number_of_lnqs2_values];
    double ff_arrays[number_of_flavors][number_of_lnz_values * number_of_lnqs2_values];

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