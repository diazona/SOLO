#include <cassert>
#include <cmath>
#include <sstream>
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
const char* GBWGluonDistribution::name() {
    return "GBW";
}


class MVIntegrationParameters {
public:
    double q; // use q instead of q^2 because that's what we need for the integrand
    double Qs2;
    size_t n; // for integrating a term of the series: this is the order of the term being integrated
    MVGluonDistribution* gdist;
    MVIntegrationParameters(MVGluonDistribution* gdist) : gdist(gdist), q(0), Qs2(0), n(0) {}
};

static double mv_gdist_integrand(double r, void* closure) {
    MVIntegrationParameters* params = (MVIntegrationParameters*)closure;
    return 0.5 * M_1_PI * r * params->gdist->S2(r*r, params->Qs2) * gsl_sf_bessel_J0(params->q * r);
}

static double mv_gdist_series_term_integrand(double r, void* closure) {
    MVIntegrationParameters* params = (MVIntegrationParameters*)closure;
    switch (params->n) {
        case 0:
            return 0.5 * M_1_PI * r * params->gdist->S2(r*r, params->Qs2);
        case 2:
            return -0.125 * M_1_PI * gsl_pow_3(r) * params->gdist->S2(r*r, params->Qs2);
        default: // a term not in the series
            return 0;
    }
}

MVGluonDistribution::MVGluonDistribution(double LambdaMV, double q2min, double q2max, double Qs2min, double Qs2max) :
 LambdaMV(LambdaMV), q2min(q2min), q2max(q2max), Qs2min(Qs2min), Qs2max(Qs2max),
 F_dist_leading_q2(NULL), F_dist_subleading_q2(NULL), F_dist(NULL),
 interp_dist_leading_q2(NULL), interp_dist_subleading_q2(NULL), interp_dist(NULL),
 q2_accel(NULL), Qs2_accel(NULL) {
    static const size_t subinterval_limit = 1000;
    double step = 1.2;
    MVIntegrationParameters params(this);
    gsl_function func;
    func.params = &params;
    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(subinterval_limit);
    
    double log_step = log(step);
    double log_q2min = log(q2min);
    double log_Qs2min = log(Qs2min);
    double log_q2max = log(q2max);
    double log_Qs2max = log(Qs2max);
    assert(log_step > 0);
    
    size_t q2_dimension = (size_t)((log_q2max - log_q2min) / log_step) + 2; // subtracting logs rather than dividing may help accuracy
    size_t Qs2_dimension = (size_t)((log_Qs2max - log_Qs2min) / log_step) + 2;
    
    while (q2_dimension < 4) { // 4 points needed for bicubic interpolation
        q2min /= step;
        log_q2min -= log_step;
        q2_dimension++;
    }
    while (Qs2_dimension < 4) { // 4 points needed for bicubic interpolation
        Qs2min /= step;
        log_Qs2min -= log_step;
        Qs2_dimension++;
    }
    
    log_q2_values = new double[q2_dimension];
    log_Qs2_values = new double[Qs2_dimension];
    
    double error; // throwaway
    
    // calculate the coefficients for the series approximation
    func.function = &mv_gdist_series_term_integrand;
    F_dist_leading_q2 = new double[Qs2_dimension]; // zeroth order term in series around q2 = 0
    F_dist_subleading_q2 = new double[Qs2_dimension]; // second order term in series around q2 = 0
    for (size_t i_Qs2 = 0; i_Qs2 < Qs2_dimension; i_Qs2++) {
        log_Qs2_values[i_Qs2] = log_Qs2min + i_Qs2 * log_step;
        params.Qs2 = exp(log_Qs2_values[i_Qs2]);
        
        params.n = 0;
        gsl_integration_qagiu(&func, 0, 0, 0.0001, subinterval_limit, workspace, F_dist_leading_q2 + i_Qs2, &error);
        
        params.n = 2;
        gsl_integration_qagiu(&func, 0, 0, 0.0001, subinterval_limit, workspace, F_dist_subleading_q2 + i_Qs2, &error);
    }
    
    // calculate the values for the 2D interpolation
    func.function = &mv_gdist_integrand;
    F_dist = new double[q2_dimension * Qs2_dimension];
    for (size_t i_q2 = 0; i_q2 < q2_dimension; i_q2++) {
        log_q2_values[i_q2] = log_q2min + i_q2 * log_step;
        params.q = exp(0.5 * log_q2_values[i_q2]);
        for (size_t i_Qs2 = 0; i_Qs2 < Qs2_dimension; i_Qs2++) {
            params.Qs2 = exp(log_Qs2_values[i_Qs2]);
            size_t index = INDEX_2D(i_q2, i_Qs2, q2_dimension, Qs2_dimension);
            gsl_integration_qagiu(&func, 0, 0, 0.0001, subinterval_limit, workspace, F_dist + index, &error);
        }
    }
    
    assert(log_q2_values[0] <= log_q2min);
    assert(log_q2_values[q2_dimension - 1] >= log_q2max);
    assert(log_Qs2_values[0] <= log_Qs2min);
    assert(log_Qs2_values[Qs2_dimension - 1] >= log_Qs2max);
    
    interp_dist_leading_q2 = gsl_interp_alloc(gsl_interp_cspline, Qs2_dimension);
    gsl_interp_init(interp_dist_leading_q2, log_Qs2_values, F_dist_leading_q2, Qs2_dimension);
    
    interp_dist_subleading_q2 = gsl_interp_alloc(gsl_interp_cspline, Qs2_dimension);
    gsl_interp_init(interp_dist_subleading_q2, log_Qs2_values, F_dist_subleading_q2, Qs2_dimension);
    
    interp_dist = interp2d_alloc(interp2d_bilinear, q2_dimension, Qs2_dimension);
    interp2d_init(interp_dist, log_q2_values, log_Qs2_values, F_dist, q2_dimension, Qs2_dimension);
    
    q2_accel = gsl_interp_accel_alloc();
    Qs2_accel = gsl_interp_accel_alloc();

    std::ostringstream s;
    s << "MV(LambdaMV = " << LambdaMV << ", q2min = " << q2min << ", q2max = " << q2max << ", Qs2min = " << Qs2min << ", Qs2max = " << Qs2max << ")";
    _name = s.str();
}

MVGluonDistribution::~MVGluonDistribution() {
    delete[] F_dist;
    F_dist = NULL;
    delete[] F_dist_leading_q2;
    F_dist_leading_q2 = NULL;
    delete[] F_dist_subleading_q2;
    F_dist_subleading_q2 = NULL;
    interp2d_free(interp_dist);
    interp_dist = NULL;
    gsl_interp_free(interp_dist_leading_q2);
    interp_dist = NULL;
    gsl_interp_free(interp_dist_subleading_q2);
    interp_dist = NULL;
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
    if (q2 > q2min) {
        return interp2d_eval(interp_dist, log_q2_values, log_Qs2_values, F_dist, log(q2), log(Qs2), q2_accel, Qs2_accel);
    }
    else {
        double c0 = gsl_interp_eval(interp_dist_leading_q2, log_Qs2_values, F_dist_leading_q2, log(Qs2), Qs2_accel);
        double c2 = gsl_interp_eval(interp_dist_subleading_q2, log_Qs2_values, F_dist_subleading_q2, log(Qs2), Qs2_accel);
        return c0 + c2 * q2;
    }
}

const char* MVGluonDistribution::name() {
    return _name.c_str();
}

std::ostream& operator<<(std::ostream& out, GluonDistribution& gdist) {
    out << gdist.name();
    return out;
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
