#include <cassert>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

static const double inf = 100;

// static double ic02;
// 
// double hypot_sq(double x, double y) {
//     return gsl_pow_2(x) + gsl_pow_2(y);
// }
// 
// double GBW_S(double r2, double Qs2) {
//     return exp(-0.25*r2*Qs2);
// }
// 
// double GBW_F(double k2, double Qs2) {
//     double iQs2 = 1.0/Qs2;
//     return M_1_PI * iQs2 * exp(-k2/Qs2);
// }
// 
// double L1_integrand(double* coords, size_t ndims, void* params) {
//     double k2 = ((double*)params)[0];
//     double Qs2 = ((double*)params)[1];
//     double k = sqrt(k2);
//     double r = coords[0];
//     double r2 = r*r;
//     double result;
//     double logk2r2 = log(k2*r2*ic02);
//     double Sr = GBW_S(r2, Qs2);
//     assert(ndims == 1);
//     
//     result = -r * gsl_sf_bessel_J0(k*r) * Sr * gsl_pow_2(logk2r2);
//     return result;
// }
// 
// double L2_integrand(double* coords, size_t ndims, void* params) {
//     double k2 = ((double*)params)[0];
//     double Qs2 = ((double*)params)[1];
//     double k = sqrt(k2);
//     double lx = coords[0];
//     double ly = coords[1];
//     double result;
//     double l2 = hypot_sq(lx, ly);
//     double logk2l2 = log(k2/l2);
//     double Fk = GBW_F(k2, Qs2);
//     double Fkl = GBW_F(hypot_sq(k - lx, ly), Qs2);
//     assert(ndims == 2);
//     
//     result = 2 * M_PI * Fk * Fkl * gsl_pow_2(logk2l2);
//     assert(gsl_finite(result));
//     return result;
// }
// 
// double L3_integrand(double* coords, size_t ndims, void* params) {
//     double k2 = ((double*)params)[0];
//     double Qs2 = ((double*)params)[1];
//     double k = sqrt(k2);
//     double lx = coords[0];
//     double ly = coords[1];
//     double lpx = coords[2];
//     double lpy = coords[3];
//     double l2 = hypot_sq(lx, ly);
//     double lp2 = hypot_sq(lpx, lpy);
//     double result;
//     assert(ndims == 4);
//     result = -4 * GBW_F(hypot_sq(k - lpx, lpy), Qs2) * GBW_F(hypot_sq(k - lx, ly), Qs2) * (lpx * lx + lpy * ly) / (l2 * lp2) * log(k2/l2);
//     assert(gsl_finite(result));
//     return result;
// }
// 
// double Lg1_integrand(double* coords, size_t ndims, void* params) {
//     double k2 = ((double*)params)[0];
//     double Qs2 = ((double*)params)[1];
//     double k = sqrt(k2);
//     double r = coords[0];
//     double r2 = r*r;
//     double result;
//     double logk2r2 = log(k2*r2*ic02);
//     double Sr = GBW_S(r2, Qs2);
//     assert(ndims == 1);
// 
//     result = -2*r*gsl_sf_bessel_J0(k*r) * gsl_pow_2(Sr) * gsl_pow_2(logk2r2);
//     assert(gsl_finite(result));
//     return result;
// }
// 
// double Lg2_integrand(double* coords, size_t ndims, void* params) {
//     double k2 = ((double*)params)[0];
//     double Qs2 = ((double*)params)[1];
//     double k = sqrt(k2);
//     double lx = coords[0];
//     double ly = coords[1];
//     double lpx = coords[2];
//     double lpy = coords[3];
//     double l2 = hypot_sq(lx, ly);
//     double lp2 = hypot_sq(lpx, lpy);
//     double result;
//     double logk2l2 = log(k2/l2);
//     assert(ndims == 4);
//     
//     result = 4*M_PI * GBW_F(hypot_sq(k - lx - lpx, -ly - lpy), Qs2) * GBW_F(hypot_sq(k - lpx, -lpy), Qs2) * GBW_F(lp2, Qs2) * gsl_pow_2(logk2l2);
//     assert(gsl_finite(result));
//     return result;
// }
// 
// double Lg3_integrand(double* coords, size_t ndims, void* params) {
//     double k2 = ((double*)params)[0];
//     double Qs2 = ((double*)params)[1];
//     double k = sqrt(k2);
//     double lx = coords[0];
//     double ly = coords[1];
//     double lpx = coords[2];
//     double lpy = coords[3];
//     double lppx = coords[4];
//     double lppy = coords[5];
//     double l2 = hypot_sq(lx, ly);
//     double lp2 = hypot_sq(lpx, lpy);
//     double lpp2 = hypot_sq(lppx, lppy);
//     double result;
//     double logk2l2 = log(k2/l2);
//     double ldotlp = (lpx * lx + lpy * ly)/(lp2*l2);
//     assert(ndims == 6);
//     
//     result = -8 * GBW_F(lpp2, Qs2) * GBW_F(hypot_sq(k - lpx - lppx, -lpy - lppy), Qs2) * GBW_F(hypot_sq(k - lx - lppx, -ly - lppy), Qs2) * ldotlp * (logk2l2);
//     assert(gsl_finite(result));
//     return result;
// }

void integrate(const unsigned int seed, double (*func)(double*, size_t, void*), const size_t ndim, const double* min, const double* max, const size_t ncalls_pre, const size_t ncalls_iter, void* params, double* result, double* error_estimate) {
    double res, err;
    const gsl_rng_type* T;
    gsl_rng *r;
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);

    gsl_monte_function function_package;
    function_package.f = func;
    function_package.dim = ndim;
    function_package.params = params;
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(ndim);
    gsl_monte_vegas_integrate (&function_package, (double*)min, (double*)max, (size_t)ndim, (size_t)ncalls_pre, r, s, &res, &err);

    do {
        gsl_monte_vegas_integrate(&function_package, (double*)min, (double*)max, (size_t)ndim, (size_t)ncalls_iter, r, s, &res, &err);
//         printf("result = % .6f sigma = % .6f chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq(s));
    } while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.1);
    
    *result = res;
    *error_estimate = err;

    gsl_monte_vegas_free (s);
    gsl_rng_free (r);
}

// integrates multiple times and averages the results
void multi_integrate(const unsigned int starting_seed, const size_t nsamples, double (*func)(double*, size_t, void*), const size_t ndim, const double* min, const double* max, const size_t ncalls_pre, const size_t ncalls_iter, void* params, double* result, double* stddev, double* error_estimate) {
    double mk = 0., mkm1 = 0., sk = 0., total_error = 0.;
    double l_result, l_error_estimate;
    size_t k = 0;
    do {
        ++k;
        integrate(starting_seed + (unsigned int)k, func, ndim, min, max, ncalls_pre, ncalls_iter, params, &l_result, &l_error_estimate);
        mkm1 = mk;
        mk += (l_result - mkm1)/(double)k;
        sk += (l_result - mkm1)*(l_result - mk);
        total_error += l_error_estimate;
    } while (k < nsamples); // doing the loop this way leaves k==nsamples after the last iteration
    *result = mk;
    *stddev = sqrt(sk)/(double)k; // standard error of the mean
    *error_estimate = total_error;
}



int main(const int argc, const char** argv) {
    const unsigned int starting_seed = 100;
    size_t nsamples = 40;
    double min0[] = {0};
    double min[] = {-inf, -inf, -inf, -inf, -inf, -inf};
    double max[] = { inf,  inf,  inf,  inf,  inf,  inf};
    double params[2]; // k2 and Qs2

    typedef double (*ftype)(double*, size_t, void*);
    ftype f[] = {L1_integrand, L2_integrand, L3_integrand, Lg1_integrand, Lg2_integrand, Lg3_integrand};
    double dims[] = {1, 2, 4, 1, 4, 6};
    double* mins[] = {min0, min, min, min0, min, min};
    double r[6];
    double s[6];
    double e[6];
    
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " filename [filename...]" << endl;
        return 2;
    }
    HardFactorParser parser;
    for (size_t i = 1; i < argc; i++) {
        try {
            parser.parse_file(argv[i], print_err);
        }
        catch (const std::exception& e) {
            cerr << e.what() << endl;
        }
        catch (const mu::ParserError& e) {
            cerr << e.GetMsg() << endl;
        }
    }
    HardFactorList hl = parser.get_hard_factors();
}