/*
 * Part of oneloopcalc
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
#include <typeinfo>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include "integrator.h"

/*
 * Here's the overall "usage map" of the code in this file:
 * 
 * The program constructs an Integrator object,
 * sets any relevant callbacks,
 * and calls integrate() on it.
 * (this happens in ResultsCalculator::calculate())
 * 
 * integrate() iterates over the IntegrationTypes, and for each, calls
 *   integrate_impl() to do the "2D" integral, which calls
 *     vegas_integrate(), which calls
 *       gsl_monte_vegas_integrate(), passing the wrapper function gsl_monte_wrapper_2D as the function to be integrated. 
 *         gsl_monte_wrapper_2D() calls
 *           Integrator::update2D() to update the values in the IntegrationContext, and then calls
 *           Integrator::evaluate_2D_integrand() to actually calculate the result
 *   Then integrate() calls
 *   integrate_impl() to do the "1D" integral, which calls
 *     vegas_integrate(), which calls
 *       gsl_monte_vegas_integrate(), passing the wrapper function gsl_monte_wrapper_1D as the function to be integrated. 
 *         gsl_monte_wrapper_1D() calls
 *           Integrator::update1D() to update the values in the IntegrationContext, and then calls
 *           Integrator::evaluate_1D_integrand() to actually calculate the result
 */

/**
 * Check whether a value is finite, i.e. not inf, -inf, or nan.
 */
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

/**
 * A wrapper function that can be passed to the GSL integration code
 * to do the 1D integrand.
 * 
 * This updates the IntegrationContext using the values in `coordinates`
 * and then evaluates the current list of hard factors. The `coordinates`
 * are interpreted as z, (xiprime if applicable), rx, ry, etc.
 */
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

/**
 * A wrapper function that can be passed to the GSL integration code
 * to do the 2D integrand.
 * 
 * This updates the IntegrationContext using the values in `coordinates`
 * and then evaluates the current list of hard factors. The `coordinates`
 * are interpreted as z, y, (xiprime if applicable), rx, ry, etc.
 */
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

/**
 * A wrapper function that calls the GSL integration routine for MISER integration.
 * 
 * @param func the integrand
 * @param dim the number of dimensions it's being integrated over
 * @param closure something to be passed to the integrand as its last argument
 * @param min the lower bounds of the integration region
 * @param max the upper bounds of the integration region
 * @param[out] p_result the result
 * @param[out] p_abserr the error bound
 * @param iterations the number of function evaluations
 * @param rng the random number generator
 * @param callback a callback to call when the integration is done
 */
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

/**
 * A wrapper function that calls the GSL integration routine for VEGAS integration.
 * 
 * @param func the integrand
 * @param dim the number of dimensions it's being integrated over
 * @param closure something to be passed to the integrand as its last argument
 * @param min the lower bounds of the integration region
 * @param max the upper bounds of the integration region
 * @param[out] p_result the result
 * @param[out] p_abserr the error bound
 * @param initial_iterations the number of function evaluations to use when first refining the grid
 * @param incremental_iterations the number of function evaluations to use in subsequent steps
 * @param rng the random number generator
 * @param callback a callback to call when the integration is done
 */
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
    switch (strategy) {
        case MC_VEGAS:
            vegas_integrate(monte_wrapper, dimensions, this, min, max, result, error, ictx.ctx->vegas_initial_iterations, ictx.ctx->vegas_incremental_iterations, rng, vegas_callback);
            break;
        case MC_MISER:
            miser_integrate(monte_wrapper, dimensions, this, min, max, result, error, ictx.ctx->miser_iterations, rng, miser_callback);
            break;
        case MC_PLAIN:
            throw "Unsupported integration method PLAIN";
        default:
            throw "Unknown integration method";
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
