#include <cassert>
#include <typeinfo>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include "integrator.h"

#define checkfinite(d) assert(gsl_finite(d))

Integrator::Integrator(Context* ctx, size_t hflen, HardFactor** hflist) : n_dipole_terms(0), n_quadrupole_terms(0), n_momentum1_terms(0), n_momentum2_terms(0), n_momentum3_terms(0), current_term_type(NONE), callback(NULL) {
    assert(hflen > 0);
    ictx = new IntegrationContext(ctx);
    // separate the hard factors provided into dipole and quadrupole terms
    for (size_t i = 0; i < hflen; i++) {
        assert(hflist[i]->get_type() != NONE);
        assert(hflist[i]->get_type() == dipole || hflist[i]->get_type() == quadrupole || hflist[i]->get_type() == momentum1 || hflist[i]->get_type() == momentum2 || hflist[i]->get_type() == momentum3);
        if (hflist[i]->get_type() == dipole) {
            n_dipole_terms++;
        }
        else if (hflist[i]->get_type() == quadrupole) {
            n_quadrupole_terms++;
        }
        else if (hflist[i]->get_type() == momentum1) {
            n_momentum1_terms++;
        }
        else if (hflist[i]->get_type() == momentum2) {
            n_momentum2_terms++;
        }
        else if (hflist[i]->get_type() == momentum3) {
            n_momentum3_terms++;
        }
    }
    assert(n_dipole_terms + n_quadrupole_terms + n_momentum1_terms + n_momentum2_terms + n_momentum3_terms == hflen);
    dipole_terms = (HardFactor**)calloc(n_dipole_terms, sizeof(HardFactor*));
    quadrupole_terms = (HardFactor**)calloc(n_quadrupole_terms, sizeof(HardFactor*));
    momentum1_terms = (HardFactor**)calloc(n_momentum1_terms, sizeof(HardFactor*));
    momentum2_terms = (HardFactor**)calloc(n_momentum2_terms, sizeof(HardFactor*));
    momentum3_terms = (HardFactor**)calloc(n_momentum3_terms, sizeof(HardFactor*));
    size_t i[6] = {0,0,0,0,0,0};
    for (; i[0] < hflen; i[0]++) {
        size_t index = i[0];
        assert(index == i[dipole] + i[quadrupole] + i[momentum1] + i[momentum2] + i[momentum3]);
        switch (hflist[index]->get_type()) {
            case dipole:
                dipole_terms[i[dipole]++] = hflist[index];
                break;
            case quadrupole:
                quadrupole_terms[i[quadrupole]++] = hflist[index];
                break;
            case momentum1:
                momentum1_terms[i[momentum1]++] = hflist[index];
                break;
            case momentum2:
                momentum2_terms[i[momentum2]++] = hflist[index];
                break;
            case momentum3:
                momentum3_terms[i[momentum3]++] = hflist[index];
                break;
            case NONE:
            case COUNT:
                assert(false);
                break;
        }
    }
    assert(i[dipole] == n_dipole_terms);
    assert(i[quadrupole] == n_quadrupole_terms);
    assert(i[momentum1] == n_momentum1_terms);
    assert(i[momentum2] == n_momentum2_terms);
    assert(i[momentum3] == n_momentum3_terms);
}

Integrator::~Integrator() {
    free(dipole_terms);
    dipole_terms = NULL;
    free(quadrupole_terms);
    quadrupole_terms = NULL;
    free(momentum1_terms);
    momentum1_terms = NULL;
    free(momentum2_terms);
    momentum2_terms = NULL;
    free(momentum3_terms);
    momentum3_terms = NULL;
    delete ictx;
}

void Integrator::update_position(double z, double y, double rx, double ry) {
    assert(current_term_type == dipole);
    ictx->update_position(z, y, rx, ry, 0, 0, 0, 0);
}

void Integrator::update_position(double z, double y, double sx, double sy, double tx, double ty) {
    assert(current_term_type == quadrupole);
    ictx->update_position(z, y, sx, sy, tx, ty, 0, 0);
}

void Integrator::update_momentum(double z, double y, double q1x, double q1y) {
    assert(current_term_type == momentum1);
    ictx->update_momentum(z, y, q1x, q1y, 0, 0, 0, 0);
}

void Integrator::update_momentum(double z, double y, double q1x, double q1y, double q2x, double q2y) {
    assert(current_term_type == momentum2);
    ictx->update_momentum(z, y, q1x, q1y, q2x, q2y, 0, 0);
}

void Integrator::update_momentum(double z, double y, double q1x, double q1y, double q2x, double q2y, double q3x, double q3y) {
    assert(current_term_type == momentum3);
    ictx->update_momentum(z, y, q1x, q1y, q2x, q2y, q3x, q3y);
}

void Integrator::evaluate_1D_integrand(double* real, double* imag) {
    double l_real = 0.0, l_imag = 0.0; // l for "local"
    double t_real, t_imag;             // t for temporary
    double log_factor = log(1 - ictx->ctx->tau / ictx->z);
    checkfinite(log_factor);
    size_t n_terms;
    HardFactor** terms;
    assert(current_term_type == dipole || current_term_type == quadrupole || current_term_type == momentum1 || current_term_type == momentum2 || current_term_type == momentum3);
    switch (current_term_type) {
        case dipole:
            n_terms = n_dipole_terms;
            terms = dipole_terms;
            break;
        case quadrupole:
            n_terms = n_quadrupole_terms;
            terms = quadrupole_terms;
            break;
        case momentum1:
            n_terms = n_momentum1_terms;
            terms = momentum1_terms;
            break;
        case momentum2:
            n_terms = n_momentum2_terms;
            terms = momentum2_terms;
            break;
        case momentum3:
            n_terms = n_momentum3_terms;
            terms = momentum3_terms;
            break;
        case NONE:
        case COUNT:
            assert(false);
            break;
    }
    assert(n_terms > 0);
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
    assert(current_term_type == dipole || current_term_type == quadrupole || current_term_type == momentum1 || current_term_type == momentum2 || current_term_type == momentum3);
    switch (current_term_type) {
        case dipole:
            n_terms = n_dipole_terms;
            terms = dipole_terms;
            break;
        case quadrupole:
            n_terms = n_quadrupole_terms;
            terms = quadrupole_terms;
            break;
        case momentum1:
            n_terms = n_momentum1_terms;
            terms = momentum1_terms;
            break;
        case momentum2:
            n_terms = n_momentum2_terms;
            terms = momentum2_terms;
            break;
        case momentum3:
            n_terms = n_momentum3_terms;
            terms = momentum3_terms;
            break;
        case NONE:
        case COUNT:
            assert(false);
            break;
    }
    assert(n_terms > 0);
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
    switch (current_term_type) {
        // set xi to 1
        case dipole:
        case quadrupole:
            ictx->update_position(ictx->z, 1, ictx->xx, ictx->xy, ictx->yx, ictx->yy, ictx->bx, ictx->by);
            break;
        case momentum1:
        case momentum2:
        case momentum3:
            ictx->update_momentum(ictx->z, 1, ictx->q1x, ictx->q1y, ictx->q2x, ictx->q2y, ictx->q3x, ictx->q3y);
            break;
        case NONE:
        case COUNT:
            assert(false);
            break;
    }
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

double gsl_monte_wrapper_1D_position(double* coordinates, size_t ncoords, void* closure) {
    double real;
    double imag;
    Integrator* integrator = (Integrator*)closure;
    assert(ncoords == 3 || ncoords == 5);
    if (ncoords == 3) {
        integrator->update_position(coordinates[0], 1, coordinates[1], coordinates[2]);
    }
    else {
        integrator->update_position(coordinates[0], 1, coordinates[1], coordinates[2], coordinates[3], coordinates[4]);
    }
    integrator->evaluate_1D_integrand(&real, &imag);
    checkfinite(real);
    checkfinite(imag);
    return real;
}

double gsl_monte_wrapper_2D_position(double* coordinates, size_t ncoords, void* closure) {
    double real;
    double imag;
    Integrator* integrator = (Integrator*)closure;
    assert(ncoords == 4 || ncoords == 6);
    if (ncoords == 4) {
        integrator->update_position(coordinates[0], coordinates[1], coordinates[2], coordinates[3]);
    }
    else {
        integrator->update_position(coordinates[0], coordinates[1], coordinates[2], coordinates[3], coordinates[4], coordinates[5]);
    }
    integrator->evaluate_2D_integrand(&real, &imag);
    checkfinite(real);
    checkfinite(imag);
    return real;
}

double gsl_monte_wrapper_1D_momentum(double* coordinates, size_t ncoords, void* closure) {
    double real;
    double imag;
    Integrator* integrator = (Integrator*)closure;
    assert(ncoords == 3 || ncoords == 5 || ncoords == 7);
    if (ncoords == 3) {
        integrator->update_momentum(coordinates[0], 1, coordinates[1], coordinates[2]);
    }
    else if (ncoords == 5) {
        integrator->update_momentum(coordinates[0], 1, coordinates[1], coordinates[2], coordinates[3], coordinates[4]);
    }
    else if (ncoords == 7) {
        integrator->update_momentum(coordinates[0], 1, coordinates[1], coordinates[2], coordinates[3], coordinates[4], coordinates[5], coordinates[6]);
    }
    integrator->evaluate_1D_integrand(&real, &imag);
    checkfinite(real);
    checkfinite(imag);
    return real;
}

double gsl_monte_wrapper_2D_momentum(double* coordinates, size_t ncoords, void* closure) {
    double real;
    double imag;
    Integrator* integrator = (Integrator*)closure;
    assert(ncoords == 4 || ncoords == 6 || ncoords == 8);
    if (ncoords == 4) {
        integrator->update_momentum(coordinates[0], coordinates[1], coordinates[2], coordinates[3]);
    }
    else if (ncoords == 6) {
        integrator->update_momentum(coordinates[0], coordinates[1], coordinates[2], coordinates[3], coordinates[4], coordinates[5]);
    }
    else if (ncoords == 8) {
        integrator->update_momentum(coordinates[0], coordinates[1], coordinates[2], coordinates[3], coordinates[4], coordinates[5], coordinates[6], coordinates[7]);
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

void Integrator::integrate_impl(double (*monte_wrapper)(double*, size_t, void*), size_t dimensions, double* min, double* max, double* result, double* error) {
    // it should already have been checked that there is at least one term of the appropriate type
    // and the type should be set appropriately
    if (integration_strategy == MC_VEGAS) {
        vegas_integrate(monte_wrapper, dimensions, this, min, max, result, error, vegas_eprint_callback);
    }
    else if (integration_strategy == MC_MISER) {
        miser_integrate(monte_wrapper, dimensions, this, min, max, result, error, miser_eprint_callback);
    }
    if (callback) {
        callback(NULL, 0, 0);
    }
}

void Integrator::integrate(double* real, double* imag, double* error) {
    double min1D[] = {ictx->ctx->tau, -inf, -inf, -inf, -inf, -inf, -inf};
    double max1D[] = {1.0, inf, inf, inf, inf, inf, inf};
    double min2D[] = {ictx->ctx->tau, ictx->ctx->tau, -inf, -inf, -inf, -inf, -inf, -inf};
    double max2D[] = {1.0, 1.0, inf, inf, inf, inf, inf, inf};
    double result = 0.0;
    double abserr = 0.0;
    double tmp_result, tmp_error;
    // cubature doesn't work because of the endpoint singularity at xi = 1
    
    // dipole
    if (n_dipole_terms > 0) {
        set_current_term_type(dipole);
        integrate_impl(gsl_monte_wrapper_2D_position, 4, min2D, max2D, &tmp_result, &tmp_error);
        result += tmp_result;
        abserr += tmp_error;
        integrate_impl(gsl_monte_wrapper_1D_position, 3, min1D, max1D, &tmp_result, &tmp_error);
        result += tmp_result;
        abserr += tmp_error;
    }
    if (n_quadrupole_terms > 0) {
        set_current_term_type(quadrupole);
        integrate_impl(gsl_monte_wrapper_2D_position, 6, min2D, max2D, &tmp_result, &tmp_error);
        result += tmp_result;
        abserr += tmp_error;
        integrate_impl(gsl_monte_wrapper_1D_position, 5, min1D, max1D, &tmp_result, &tmp_error);
        result += tmp_result;
        abserr += tmp_error;
    }
    if (n_momentum1_terms > 0) {
        set_current_term_type(momentum1);
        integrate_impl(gsl_monte_wrapper_2D_momentum, 4, min2D, max2D, &tmp_result, &tmp_error);
        result += tmp_result;
        abserr += tmp_error;
        integrate_impl(gsl_monte_wrapper_1D_momentum, 3, min1D, max1D, &tmp_result, &tmp_error);
        result += tmp_result;
        abserr += tmp_error;
    }
    if (n_momentum2_terms > 0) {
        set_current_term_type(momentum2);
        integrate_impl(gsl_monte_wrapper_2D_momentum, 6, min2D, max2D, &tmp_result, &tmp_error);
        result += tmp_result;
        abserr += tmp_error;
        integrate_impl(gsl_monte_wrapper_1D_momentum, 5, min1D, max1D, &tmp_result, &tmp_error);
        result += tmp_result;
        abserr += tmp_error;
    }
    if (n_momentum3_terms > 0) {
        set_current_term_type(momentum3);
        integrate_impl(gsl_monte_wrapper_2D_momentum, 8, min2D, max2D, &tmp_result, &tmp_error);
        result += tmp_result;
        abserr += tmp_error;
        integrate_impl(gsl_monte_wrapper_1D_momentum, 7, min1D, max1D, &tmp_result, &tmp_error);
        result += tmp_result;
        abserr += tmp_error;
    }

    *real = result;
    *imag = 0;
    *error = abserr;
}

