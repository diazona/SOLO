#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include "integrator.h"

#define checkfinite(d) assert(gsl_finite(d))

Integrator::Integrator(Context* ctx, size_t hflen, HardFactor** hflist) : n_dipole_terms(0), n_quadrupole_terms(0), current_term_type(NONE), callback(NULL) {
    size_t i, j;
    ictx = new IntegrationContext(ctx);
    // separate the hard factors provided into dipole and quadrupole terms
    for (i = 0; i < hflen; i++) {
        assert(hflist[i]->get_type() != NONE);
        assert(hflist[i]->get_type() == dipole || hflist[i]->get_type() == quadrupole);
        if (hflist[i]->get_type() == dipole) {
            n_dipole_terms++;
        }
        else if (hflist[i]->get_type() == quadrupole) {
            n_quadrupole_terms++;
        }
    }
    assert(n_dipole_terms + n_quadrupole_terms == hflen);
    dipole_terms = (HardFactor**)calloc(n_dipole_terms, sizeof(HardFactor*));
    quadrupole_terms = (HardFactor**)calloc(n_quadrupole_terms, sizeof(HardFactor*));
    for (i = 0, j = 0; i + j < hflen;) {
        if (hflist[i + j]->get_type() == dipole) {
            dipole_terms[i] = hflist[i + j];
            i++;
        }
        else {
            assert(hflist[i + j]->get_type() == quadrupole);
            quadrupole_terms[j] = hflist[i + j];
            j++;
        }
    }
    assert(i == n_dipole_terms);
    assert(j == n_quadrupole_terms);
}

Integrator::~Integrator() {
    free(dipole_terms);
    dipole_terms = NULL;
    free(quadrupole_terms);
    quadrupole_terms = NULL;
    delete ictx;
}


void Integrator::evaluate_1D_integrand(double* real, double* imag) {
    double l_real = 0.0, l_imag = 0.0; // l for "local"
    double t_real, t_imag;             // t for temporary
    double log_factor = log(1 - ictx->ctx->tau / ictx->z);
    checkfinite(log_factor);
    size_t n_terms;
    HardFactor** terms;
    assert(current_term_type == dipole || current_term_type == quadrupole);
    if (current_term_type == dipole) {
        n_terms = n_dipole_terms;
        terms = dipole_terms;
    }
    else {
        n_terms = n_quadrupole_terms;
        terms = quadrupole_terms;
    }
    if (n_terms == 0) {
        *real = 0;
        *imag = 0;
        return;
    }
    assert(ictx->xi == 1.0d);
    for (size_t i = 0; i < n_terms; i++) {
        terms[i]->Fs(ictx, &t_real, &t_imag);
        checkfinite(t_real);
        checkfinite(t_imag);
        l_real += t_real * log_factor;
        l_imag += t_imag * log_factor;
        terms[i]->Fd(ictx, &t_real, &t_imag);
        checkfinite(t_real);
        checkfinite(t_imag);
        l_real += t_real;
        l_imag += t_imag;
    }
    if (callback) {
        callback(ictx, l_real, l_imag);
    }
    *real = l_real;
    *imag = l_imag;
}

void Integrator::evaluate_2D_integrand(double* real, double* imag) {
    double l_real = 0.0, l_imag = 0.0; // l for "local"
    double t_real, t_imag;             // t for temporary
    double jacobian =  (1 - ictx->ctx->tau / ictx->z) / (1 - ictx->ctx->tau); // Jacobian from y to xi
    checkfinite(jacobian);
    double xi_factor = 1.0 / (1 - ictx->xi);
    checkfinite(xi_factor);
    size_t n_terms;
    HardFactor** terms;
    assert(current_term_type == dipole || current_term_type == quadrupole);
    if (current_term_type == dipole) {
        n_terms = n_dipole_terms;
        terms = dipole_terms;
    }
    else {
        n_terms = n_quadrupole_terms;
        terms = quadrupole_terms;
    }
    if (n_terms == 0) {
        *real = 0;
        *imag = 0;
        return;
    }
    for (size_t i = 0; i < n_terms; i++) {
        terms[i]->Fs(ictx, &t_real, &t_imag);
        checkfinite(t_real);
        checkfinite(t_imag);
        l_real += t_real * xi_factor;
        l_imag += t_imag * xi_factor;
        terms[i]->Fn(ictx, &t_real, &t_imag);
        checkfinite(t_real);
        checkfinite(t_imag);
        l_real += t_real;
        l_imag += t_imag;
    }
    if (callback) {
        callback(ictx, jacobian * l_real, jacobian * l_imag);
    }
    ictx->update(ictx->z, 1, ictx->xx, ictx->xy, ictx->yx, ictx->yy, ictx->bx, ictx->by);
    for (size_t i = 0; i < n_terms; i++) {
        terms[i]->Fs(ictx, &t_real, &t_imag);
        checkfinite(t_real);
        checkfinite(t_imag);
        l_real -= t_real * xi_factor;
        l_imag -= t_imag * xi_factor;
    }
    if (callback) {
        callback(ictx, jacobian * l_real, jacobian * l_imag);
    }
    *real = jacobian * l_real;
    *imag = jacobian * l_imag;
}

double gsl_monte_wrapper_1D(double* coordinates, size_t ncoords, void* closure) {
    double real;
    double imag;
    Integrator* integrator = (Integrator*)closure;
    assert(ncoords == 3 || ncoords == 5);
    if (ncoords == 3) {
        integrator->update(coordinates[0], 1, coordinates[1], coordinates[2]);
    }
    else {
        integrator->update(coordinates[0], 1, coordinates[1], coordinates[2], coordinates[3], coordinates[4]);
    }
    integrator->evaluate_1D_integrand(&real, &imag);
    checkfinite(real);
    checkfinite(imag);
    return real;
}

double gsl_monte_wrapper_2D(double* coordinates, size_t ncoords, void* closure) {
    double real;
    double imag;
    Integrator* integrator = (Integrator*)closure;
    assert(ncoords == 4 || ncoords == 6);
    if (ncoords == 4) {
        integrator->update(coordinates[0], coordinates[1], coordinates[2], coordinates[3]);
    }
    else {
        integrator->update(coordinates[0], coordinates[1], coordinates[2], coordinates[3], coordinates[4], coordinates[5]);
    }
    integrator->evaluate_2D_integrand(&real, &imag);
    checkfinite(real);
    checkfinite(imag);
    return real;
}

void miser_integrate(double (*func)(double*, size_t, void*), size_t dim, void* closure, double* min, double* max, double* p_result, double* p_abserr,
               void (*callback)(double*, double*, gsl_monte_miser_state*)) {
    gsl_monte_function f;
    f.f = func;
    f.dim = dim;
    f.params = closure;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
    gsl_monte_miser_state* s = gsl_monte_miser_alloc(dim);
    gsl_monte_miser_integrate(&f, min, max, dim, 10000000, rng, s, p_result, p_abserr);
    checkfinite(*p_result);
    checkfinite(*p_abserr);
    if (callback) {
        (*callback)(p_result, p_abserr, s);
    }
    gsl_monte_miser_free(s);
    s = NULL;
    gsl_rng_free(rng);
    rng = NULL;
}

void vegas_integrate(double (*func)(double*, size_t, void*), size_t dim, void* closure, double* min, double* max, double* p_result, double* p_abserr,
               void (*callback)(double*, double*, gsl_monte_vegas_state*)) {
    gsl_monte_function f;
    f.f = func;
    f.dim = dim;
    f.params = closure;

    gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
    gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(dim);
    gsl_monte_vegas_integrate(&f, min, max, dim, 10000, rng, s, p_result, p_abserr);
    checkfinite(*p_result);
    checkfinite(*p_abserr);
    if (callback) {
        (*callback)(p_result, p_abserr, s);
    }
    do {
        gsl_monte_vegas_integrate(&f, min, max, dim, 100000, rng, s, p_result, p_abserr);
        checkfinite(*p_result);
        checkfinite(*p_abserr);
        if (callback) {
            (*callback)(p_result, p_abserr, s);
        }
    } while (*p_abserr > 0 && fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.2);
    gsl_monte_vegas_free(s);
    s = NULL;
    gsl_rng_free(rng);
    rng = NULL;
}

void vegas_eprint_callback(double* p_result, double* p_abserr, gsl_monte_vegas_state* s) {
    cerr << "VEGAS output: " << *p_result << " err: " << *p_abserr << " chisq:" << gsl_monte_vegas_chisq(s) << endl;
}
void miser_eprint_callback(double* p_result, double* p_abserr, gsl_monte_miser_state* s) {
    cerr << "MISER output: " << *p_result << " err: " << *p_abserr << endl;
}

static double inf = 10;

void Integrator::integrate(double* real, double* imag, double* error) {
    double tmp_result;
    double tmp_error;
    double min1D[] = {ictx->ctx->tau, -inf, -inf, -inf, -inf};
    double max1D[] = {1.0, inf, inf, inf, inf};
    double min2D[] = {ictx->ctx->tau, ictx->ctx->tau, -inf, -inf, -inf, -inf};
    double max2D[] = {1.0, 1.0, inf, inf, inf, inf};
    double result = 0.0;
    double absvar = 0.0;
    // cubature doesn't work because of the endpoint singularity at xi = 1
    
    // dipole
    if (n_dipole_terms > 0) {
        set_current_term_type(dipole);
        // 2D integral
        if (integration_strategy == MC_VEGAS) {
            vegas_integrate(gsl_monte_wrapper_2D, 4, this, min2D, max2D, &tmp_result, &tmp_error, vegas_eprint_callback);
        }
        else if (integration_strategy == MC_MISER) {
            miser_integrate(gsl_monte_wrapper_2D, 4, this, min2D, max2D, &tmp_result, &tmp_error, miser_eprint_callback);
        }
        if (callback) {
            callback(NULL, 0, 0);
        }
        result += tmp_result;
        absvar += tmp_error * tmp_error;
        // 1D integral
        if (integration_strategy == MC_VEGAS) {
            vegas_integrate(gsl_monte_wrapper_1D, 3, this, min1D, max1D, &tmp_result, &tmp_error, vegas_eprint_callback);
        }
        else if (integration_strategy == MC_MISER) {
            miser_integrate(gsl_monte_wrapper_1D, 3, this, min1D, max1D, &tmp_result, &tmp_error, miser_eprint_callback);
        }
        if (callback) {
            callback(NULL, 0, 0);
        }
        result += tmp_result;
        absvar += tmp_error * tmp_error;
    }
    if (n_quadrupole_terms > 0) {
        // quadrupole
        set_current_term_type(quadrupole);
        // 2D integral
        if (integration_strategy == MC_VEGAS) {
            vegas_integrate(gsl_monte_wrapper_2D, 6, this, min2D, max2D, &tmp_result, &tmp_error, vegas_eprint_callback);
        }
        else if (integration_strategy == MC_MISER) {
            miser_integrate(gsl_monte_wrapper_2D, 6, this, min2D, max2D, &tmp_result, &tmp_error, miser_eprint_callback);
        }
        if (callback) {
            callback(NULL, 0, 0);
        }
        result += tmp_result;
        absvar += tmp_error * tmp_error;
        // 1D integral
        if (integration_strategy == MC_VEGAS) {
            vegas_integrate(gsl_monte_wrapper_1D, 5, this, min1D, max1D, &tmp_result, &tmp_error, vegas_eprint_callback);
        }
        else if (integration_strategy == MC_MISER) {
            miser_integrate(gsl_monte_wrapper_1D, 5, this, min1D, max1D, &tmp_result, &tmp_error, miser_eprint_callback);
        }
        if (callback) {
            callback(NULL, 0, 0);
        }
        result += tmp_result;
        absvar += tmp_error * tmp_error;
    }

    *real = result;
    *imag = 0;
    *error = sqrt(absvar);
}

