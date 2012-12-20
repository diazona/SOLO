#include <cassert>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include "gluondist.h"
#include "interp2d.h"

double GBWGluonDistribution::S2(double r2, double Qs2) {
    return exp(-0.25 * r2 * Qs2);
}
double GBWGluonDistribution::S4(double r2, double s2, double t2, double Qs2) {
    return exp(-0.25 * Qs2 * (s2 + t2));
}
double GBWGluonDistribution::F(double q2, double Qs2) {
    return M_1_PI * exp(-q2/Qs2) / Qs2;
}

class MVIntegrationParameters {
public:
    double q; // use q instead of q^2 because that's what we need for the integration
    double Qs2;
    MVGluonDistribution* gdist;
    MVIntegrationParameters(MVGluonDistribution* gdist) : gdist(gdist), q(0), Qs2(0) {}
};

static double mv_gdist_integrand(double r, void* closure) {
    MVIntegrationParameters* params = (MVIntegrationParameters*)closure;
    return 0.5 * M_1_PI * r * params->gdist->S2(r*r, params->Qs2) * gsl_sf_bessel_J0(params->q * r);
}

static const size_t step = 2;

MVGluonDistribution::MVGluonDistribution(double LambdaMV, double q2min, double q2max, double Qs2min, double Qs2max) :
 LambdaMV(LambdaMV), q2min(q2min), q2max(q2max), Qs2min(Qs2min), Qs2max(Qs2max), F_dist(NULL), F_error(NULL), interp_dist(NULL), interp_error(NULL), q2_accel(NULL), Qs2_accel(NULL) {
    static const size_t subinterval_limit = 1000;
    MVIntegrationParameters params(this);
    gsl_function func;
    func.function = &mv_gdist_integrand;
    func.params = &params;
    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(subinterval_limit);
    
    const double log_step = log(step);
    const double log_q2min = log(q2min);
    const double log_Qs2min = log(Qs2min);
    size_t q2_dimension = (size_t)((log(q2max) - log_q2min) / log_step) + 1; // subtracting logs rather than dividing may help accuracy
    size_t Qs2_dimension = (size_t)((log(Qs2max) - log_Qs2min) / log_step) + 1;
    
    log_q2_values = new double[q2_dimension];
    log_Qs2_values = new double[Qs2_dimension];
    
    F_dist = new double[q2_dimension * Qs2_dimension];
    F_error = new double[q2_dimension * Qs2_dimension];
    for (size_t i_q2 = 0; i_q2 < q2_dimension; i_q2++) {
        log_q2_values[i_q2] = log_q2min + i_q2 * log_step;
        params.q = exp(0.5 * log_q2_values[i_q2]);
        for (size_t i_Qs2 = 0; i_Qs2 < Qs2_dimension; i_Qs2++) {
            if (i_q2 == 0) { // only have to set the values in the array once
                log_Qs2_values[i_Qs2] = log_Qs2min + i_Qs2 * log_step;
            }
            params.Qs2 = exp(log_Qs2_values[i_Qs2]);
            size_t index = INDEX_2D(i_q2, i_Qs2, q2_dimension, Qs2_dimension);
            gsl_integration_qagiu(&func, 0, 0, 0.0001, subinterval_limit, workspace, F_dist + index, F_error + index);
        }
    }
    
    interp_dist = interp2d_alloc(interp2d_bilinear, q2_dimension, Qs2_dimension);
    interp2d_init(interp_dist, log_q2_values, log_Qs2_values, F_dist, q2_dimension, Qs2_dimension);
    
    interp_error = interp2d_alloc(interp2d_bilinear, q2_dimension, Qs2_dimension);
    interp2d_init(interp_error, log_q2_values, log_Qs2_values, F_error, q2_dimension, Qs2_dimension);

    q2_accel = gsl_interp_accel_alloc();
    Qs2_accel = gsl_interp_accel_alloc();
}

MVGluonDistribution::~MVGluonDistribution() {
    delete[] F_dist;
    F_dist = NULL;
    delete[] F_error;
    F_error = NULL;
    interp2d_free(interp_dist);
    interp_dist = NULL;
    interp2d_free(interp_error);
    interp_error = NULL;
    gsl_interp_accel_free(q2_accel);
    q2_accel = NULL;
    gsl_interp_accel_free(Qs2_accel);
    Qs2_accel = NULL;
}

double MVGluonDistribution::S2(double r2, double Qs2) {
    return pow(M_E + 1.0 / (sqrt(r2) * LambdaMV), -0.25 * r2 * Qs2);
}
double MVGluonDistribution::S4(double r2, double s2, double t2, double Qs2) {
    return S2(s2, Qs2) * S2(t2, Qs2);
}
double MVGluonDistribution::F(double q2, double Qs2) {
    return interp2d_eval(interp_dist, log_q2_values, log_Qs2_values, F_dist, log(q2), log(Qs2), q2_accel, Qs2_accel);
}

#ifdef GLUON_DIST_DRIVER
#include <cstdlib>
#include <iostream>

void MVGluonDistribution::write_grid() {
    size_t q2_dimension = interp_dist->xsize; // shouldn't really access these properties directly
    size_t Qs2_dimension = interp_dist->ysize;
    std::cout << "q2\tQs2\tF" << std::endl;
    for (size_t i_q2 = 0; i_q2 < q2_dimension; i_q2++) {
        for (size_t i_Qs2 = 0; i_Qs2 < Qs2_dimension; i_Qs2++) {
            std::cout << exp(log_q2_values[i_q2]) << "\t"
                      << exp(log_Qs2_values[i_Qs2]) << "\t"
                      << F_dist[INDEX_2D(i_q2, i_Qs2, q2_dimension, Qs2_dimension)] << std::endl;
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 6) {
        std::cerr << "Needs 5 arguments: LambdaMV, q2min, q2max, Qs2min, Qs2max" << std::endl;
        exit(1);
    }
    double LambdaMV = strtod(argv[1], NULL);
    double q2min = strtod(argv[2], NULL);
    double q2max = strtod(argv[3], NULL);
    double Qs2min = strtod(argv[4], NULL);
    double Qs2max = strtod(argv[5], NULL);
    MVGluonDistribution gdist(LambdaMV, q2min, q2max, Qs2min, Qs2max);
    std::cin.peek();
    if (std::cin.eof()) {
        // write out the grid
        gdist.write_grid();
    }
    else {
        double q2, Qs2;
        while (!std::cin.eof()) {
            std::cin >> q2 >> Qs2;
            std::cout << q2 << "\t"<< Qs2 << "\t" << gdist.F(q2, Qs2) << std::endl;
        }
    }
}
#endif
