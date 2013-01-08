#include <cassert>
#include <typeinfo>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include "integrator.h"

#define checkfinite(d) assert(gsl_finite(d))

Integrator::Integrator(const Context* ctx, const ThreadLocalContext* tlctx, const integration_strategy strategy, const HardFactorList hflist) :
  ictx(ctx, tlctx), strategy(strategy), current_type(DipoleIntegrationType::get_instance()),
  callback(NULL), miser_callback(NULL), vegas_callback(NULL) {
    assert(hflist.size() > 0);
#ifndef NDEBUG
    size_t total1 = 0;
#endif
    // separate the hard factors provided into dipole and quadrupole terms
    for (HardFactorList::const_iterator it = hflist.begin(); it != hflist.end(); it++) {
        const HardFactorTerm* const* l_terms = (*it)->get_terms();
        for (size_t i = 0; i < (*it)->get_term_count(); i++) {
            const IntegrationType* type = l_terms[i]->get_type();
            terms[type].push_back(l_terms[i]);
#ifndef NDEBUG
            total1++;
#endif
        }
    }
#ifndef NDEBUG
    size_t total2 = 0;
    for (HardFactorTypeMap::iterator it = terms.begin(); it != terms.end(); it++) {
        total2 += (*it).second.size();
    }
    assert(total1 == total2);
#endif
}

Integrator::~Integrator() {
}

void Integrator::update1D(const double* values) {
    current_type->update(ictx, 1, values);
}

void Integrator::update2D(const double* values) {
    current_type->update(ictx, 2, values);
}

void Integrator::evaluate_1D_integrand(double* real, double* imag) {
    double l_real = 0.0, l_imag = 0.0; // l for "local"
    double t_real, t_imag;             // t for temporary
    double log_factor = log(1 - ictx.ctx->tau / ictx.z);
    checkfinite(log_factor);
    HardFactorTermList& current_terms = terms[current_type];
    assert(current_terms.size() > 0);
    assert(ictx.xi == 1.0d);
    for (HardFactorTermList::const_iterator it = current_terms.begin(); it != current_terms.end(); it++) {
        const HardFactorTerm* h = (*it);
        h->Fs(&ictx, &t_real, &t_imag);
        checkfinite(t_real);
        checkfinite(t_imag);
        l_real += t_real * log_factor;
        l_imag += t_imag * log_factor;
        h->Fd(&ictx, &t_real, &t_imag);
        checkfinite(t_real);
        checkfinite(t_imag);
        l_real += t_real;
        l_imag += t_imag;
    }
    if (callback) {
        callback(&ictx, l_real, l_imag);
    }
    *real = l_real;
    *imag = l_imag;
}

void Integrator::evaluate_2D_integrand(double* real, double* imag) {
    double l_real = 0.0, l_imag = 0.0; // l for "local"
    double t_real, t_imag;             // t for temporary
    double jacobian =  (1 - ictx.ctx->tau / ictx.z) / (1 - ictx.ctx->tau); // Jacobian from y to xi
    checkfinite(jacobian);
    double xi_factor = 1.0 / (1 - ictx.xi);
    checkfinite(xi_factor);
    HardFactorTermList& current_terms = terms[current_type];
    assert(current_terms.size() > 0);
    for (HardFactorTermList::const_iterator it = current_terms.begin(); it != current_terms.end(); it++) {
        const HardFactorTerm* h = (*it);
        h->Fs(&ictx, &t_real, &t_imag);
        checkfinite(t_real);
        checkfinite(t_imag);
        l_real += t_real * xi_factor;
        l_imag += t_imag * xi_factor;
        h->Fn(&ictx, &t_real, &t_imag);
        checkfinite(t_real);
        checkfinite(t_imag);
        l_real += t_real;
        l_imag += t_imag;
    }
    if (callback) {
        callback(&ictx, jacobian * l_real, jacobian * l_imag);
    }
    ictx.set_xi_to_1();
    for (HardFactorTermList::const_iterator it = current_terms.begin(); it != current_terms.end(); it++) {
        const HardFactorTerm* h = (*it);
        h->Fs(&ictx, &t_real, &t_imag);
        checkfinite(t_real);
        checkfinite(t_imag);
        l_real -= t_real * xi_factor;
        l_imag -= t_imag * xi_factor;
    }
    if (callback) {
        callback(&ictx, jacobian * l_real, jacobian * l_imag);
    }
    *real = jacobian * l_real;
    *imag = jacobian * l_imag;
}

double gsl_monte_wrapper_1D(double* coordinates, size_t ncoords, void* closure) {
    double real;
    double imag;
    Integrator* integrator = (Integrator*)closure;
    integrator->update1D(coordinates);
    integrator->evaluate_1D_integrand(&real, &imag);
    checkfinite(real);
    checkfinite(imag);
    return real;
}

double gsl_monte_wrapper_2D(double* coordinates, size_t ncoords, void* closure) {
    double real;
    double imag;
    Integrator* integrator = (Integrator*)closure;
    integrator->update2D(coordinates);
    integrator->evaluate_2D_integrand(&real, &imag);
    checkfinite(real);
    checkfinite(imag);
    return real;
}

void miser_integrate(double (*func)(double*, size_t, void*), size_t dim, void* closure, double* min, double* max, double* p_result, double* p_abserr,
               size_t iterations, gsl_rng* rng, void (*callback)(double*, double*, gsl_monte_miser_state*)) {
    gsl_monte_function f;
    f.f = func;
    f.dim = dim;
    f.params = closure;

    gsl_monte_miser_state* s = gsl_monte_miser_alloc(dim);
    gsl_monte_miser_integrate(&f, min, max, dim, iterations, rng, s, p_result, p_abserr);
    checkfinite(*p_result);
    checkfinite(*p_abserr);
    if (callback) {
        (*callback)(p_result, p_abserr, s);
    }
    gsl_monte_miser_free(s);
    s = NULL;
}

void vegas_integrate(double (*func)(double*, size_t, void*), size_t dim, void* closure, double* min, double* max, double* p_result, double* p_abserr,
               size_t initial_iterations, size_t incremental_iterations, gsl_rng* rng, void (*callback)(double*, double*, gsl_monte_vegas_state*)) {
    gsl_monte_function f;
    f.f = func;
    f.dim = dim;
    f.params = closure;

    gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(dim);
    gsl_monte_vegas_integrate(&f, min, max, dim, initial_iterations, rng, s, p_result, p_abserr);
    checkfinite(*p_result);
    checkfinite(*p_abserr);
    if (callback) {
        (*callback)(p_result, p_abserr, s);
    }
    do {
        gsl_monte_vegas_integrate(&f, min, max, dim, incremental_iterations, rng, s, p_result, p_abserr);
        checkfinite(*p_result);
        checkfinite(*p_abserr);
        if (callback) {
            (*callback)(p_result, p_abserr, s);
        }
    } while (*p_abserr > 0 && fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.2);
    gsl_monte_vegas_free(s);
    s = NULL;
}

void Integrator::integrate_impl(size_t core_dimensions, double* result, double* error) {
    // it should already have been checked that there is at least one term of the appropriate type
    // and the type should be set appropriately
    assert(core_dimensions == 1 || core_dimensions == 2);
    size_t dimensions = core_dimensions + current_type->extra_dimensions;
    double (*monte_wrapper)(double*, size_t, void*);
    double min[10];
    double max[10];
    if (core_dimensions == 2) {
        monte_wrapper = gsl_monte_wrapper_2D;
    }
    else {
        monte_wrapper = gsl_monte_wrapper_1D;
    }
    current_type->fill_min(ictx, core_dimensions, min);
    current_type->fill_max(ictx, core_dimensions, max);
    
    gsl_rng* rng = gsl_rng_alloc(ictx.ctx->pseudorandom_generator_type);
    gsl_rng_set(rng, ictx.ctx->pseudorandom_generator_seed);
    if (strategy == MC_VEGAS) {
        vegas_integrate(monte_wrapper, dimensions, this, min, max, result, error, ictx.ctx->vegas_initial_iterations, ictx.ctx->vegas_incremental_iterations, rng, vegas_callback);
    }
    else if (strategy == MC_MISER) {
        miser_integrate(monte_wrapper, dimensions, this, min, max, result, error, ictx.ctx->miser_iterations, rng, miser_callback);
    }
    gsl_rng_free(rng);
    rng = NULL;
    if (callback) {
        callback(NULL, 0, 0);
    }
}

void Integrator::integrate(double* real, double* imag, double* error) {
    double result = 0.0;
    double abserr = 0.0;
    double tmp_result, tmp_error;
    // cubature doesn't work because of the endpoint singularity at xi = 1
    
    // dipole
    for (HardFactorTypeMap::iterator it = terms.begin(); it != terms.end(); it++) {
        assert(it->second.size() > 0);
        set_current_integration_type(it->first);
        integrate_impl(2, &tmp_result, &tmp_error);
        result += tmp_result;
        abserr += tmp_error;
        integrate_impl(1, &tmp_result, &tmp_error);
        result += tmp_result;
        abserr += tmp_error;
    }

    *real = result;
    *imag = 0;
    *error = abserr;
}
