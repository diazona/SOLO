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
#include <gsl/gsl_qrng.h>
#include "integrator.h"
#include "quasimontecarlo.h"
#include "cubature.h"

/*
 * Here's the overall "usage map" of the code in this file:
 *
 * The program constructs an Integrator object,
 * sets any relevant callbacks,
 * and calls integrate() on it.
 * (this happens in ResultsCalculator::calculate())
 *
 * integrate() iterates over the IntegrationRegions, and for each, calls
 *   integrate_impl() to do the "2D" integral, which calls
 *     vegas_integrate(), which calls
 *       gsl_monte_vegas_integrate(), passing the wrapper function gsl_monte_wrapper as the function to be integrated.
 *         gsl_monte_wrapper() calls
 *           IntegrationContext::update() to update the values in the IntegrationContext, and then calls
 *           Integrator::evaluate_integrand() to actually calculate the result
 *   Then integrate() calls
 *   integrate_impl() to do the "1D" integral, which calls
 *     vegas_integrate(), which calls
 *       gsl_monte_vegas_integrate(), passing the wrapper function gsl_monte_wrapper as the function to be integrated.
 *         gsl_monte_wrapper() calls
 *           IntegrationContext::update() to update the values in the IntegrationContext, and then calls
 *           Integrator::evaluate_integrand() to actually calculate the result
 */

/**
 * Check whether a value is finite, i.e. not inf, -inf, or nan.
 */
#define checkfinite(d) assert(gsl_finite(d))

bool compare_integration_region_pointers(const IntegrationRegion* a, const IntegrationRegion* b) {
    return a < b;
}

Integrator::Integrator(
    const Context& ctx,
    const ThreadLocalContext& tlctx,
    const HardFactorList& hflist,
    const double xg_min,
    const double xg_max) :
  ictx(ctx, tlctx),
  current_integration_region(NULL),
  xi_preintegrated_term(false),
  terms(compare_integration_region_pointers),
  xg_min(xg_min),
  xg_max(xg_max),
  callback(NULL),
  cubature_callback(NULL),
  miser_callback(NULL),
  vegas_callback(NULL),
  quasi_callback(NULL) {
    assert(hflist.size() > 0);
#ifndef NDEBUG
    size_t total1 = 0;
#endif
    // separate the hard factors provided into dipole and quadrupole etc. terms
    for (HardFactorList::const_iterator it = hflist.begin(); it != hflist.end(); it++) {
        const HardFactor* p_hf = *it;
        const HardFactorTerm* const* l_terms = p_hf->get_terms();
        for (size_t i = 0; i < p_hf->get_term_count(); i++) {
            const HardFactorTerm* term = l_terms[i];
            if (ictx.ctx.exact_kinematics && term->get_order() == HardFactor::MIXED) {
                throw KinematicSchemeMismatchException(*term);
            }
            const IntegrationRegion* region = term->get_integration();
            terms[region].push_back(term);
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

static inline bool xg_in_range(const double xg, const double xg_min, const double xg_max) {
    assert(xg_min <= xg_max);
    return xg > xg_min && xg <= xg_max;
}

void Integrator::evaluate_integrand(double* real, double* imag) {
    if (!xg_in_range(ictx.xa, xg_min, xg_max)) {
        *real = *imag = 0.0;
        return;
    }
    double l_real = 0.0, l_imag = 0.0; // l for "local"
    double t_real, t_imag;             // t for temporary
    HardFactorTermList& current_terms = terms[current_integration_region];
    assert(current_terms.size() > 0);
    if (xi_preintegrated_term) {
        // This evaluates the [Fs(1) ln(1 - tau/z) + Fd(1)] term
        assert(ictx.xi == 1.0);
        double log_factor = log(1 - ictx.ctx.tau / ictx.z);
        checkfinite(log_factor);
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
    }
    else {
        // This branch evaluates the [Fs(xi) - Fs(1)] / (1 - xi) + Fn(xi) terms
        double xi_factor = 1.0 / (1 - ictx.xi);
        checkfinite(xi_factor);
        for (HardFactorTermList::const_iterator it = current_terms.begin(); it != current_terms.end(); it++) {
            const HardFactorTerm* h = (*it);
            if (ictx.ctx.exact_kinematics) {
                // double check that there are no mixed-order hard factors when using exact kinematics
                assert(h->get_order() == HardFactor::LO || h->get_order() == HardFactor::NLO);
            }

            h->Fs(&ictx, &t_real, &t_imag);
            checkfinite(t_real);
            checkfinite(t_imag);
            if (h->get_order() == HardFactor::LO) {
                /* Leading order hard factors are supposed to only have Fd, not Fs or Fn.
                * If this assumption is violated, it could break things.
                *
                * Exercise for the reader: what things? (mwahaha)
                */
                assert(t_real == 0);
                assert(t_imag == 0);
            }
            l_real += t_real * xi_factor;
            l_imag += t_imag * xi_factor;

            h->Fn(&ictx, &t_real, &t_imag);
            checkfinite(t_real);
            checkfinite(t_imag);
            if (h->get_order() == HardFactor::LO) {
                assert(t_real == 0);
                assert(t_imag == 0);
            }
            l_real += t_real;
            l_imag += t_imag;
        }
        ictx.xi = 1;
        /* TODO put the same thing here used below in cubature_wrapper
         */
        for (HardFactorTermList::const_iterator it = current_terms.begin(); it != current_terms.end(); it++) {
            const HardFactorTerm* h = (*it);
            h->Fs(&ictx, &t_real, &t_imag);
            checkfinite(t_real);
            checkfinite(t_imag);
            l_real -= t_real * xi_factor;
            l_imag -= t_imag * xi_factor;
        }
    }
    *real = l_real;
    *imag = l_imag;
    if (callback) {
        callback(&ictx, *real, *imag);
    }
    checkfinite(*real);
    checkfinite(*imag);
}

/**
 * A wrapper function that can be passed to the cubature integration code.
 *
 * This updates the IntegrationContext using the values in `coordinates`
 * and then evaluates the current list of hard factors. The `coordinates`
 * are interpreted as z, (y if applicable), rx, ry, etc.
 */
void cubature_wrapper(unsigned int ncoords, const double* coordinates, void* closure, unsigned int nresults, double* results) {
    double real, imag, jacobian;
    Integrator* integrator = static_cast<Integrator*>(closure);
    integrator->current_integration_region->update(integrator->ictx, coordinates);
    /* TODO put here something that implements the following pseudocode:
     *
     * if (current_integration_region->position_like) {
     *     ictx.recalculate_everything_from_position(current_integration_region->is_quadrupole, current_integration_region->divide_xi);
     * }
     * else {
     *     ictx.recalculate_everything_from_momentum(current_integration_region->momentum_dimensions, current_integration_region->exact, current_integration_region->divide_xi);
     * }
     *
     * NOTE: an idea is to add a field evaluation_profile to Integrator, initialized from a profile=... line
     * in the hard factor definitions, and use that to store the settings for exact and divide_xi
     */
    // computing the Jacobian here allows the method to access the untransformed coordinates
    jacobian = integrator->current_integration_region->jacobian(integrator->ictx);
    integrator->evaluate_integrand(&real, &imag);
    assert(nresults == 1 || nresults == 2);
    results[0] = real * jacobian;
    if (nresults == 2) {
        results[1] = imag * jacobian;
    }
}

/**
 * A wrapper function that can be passed to the GSL integration code.
 *
 * This updates the IntegrationContext using the values in `coordinates`
 * and then evaluates the current list of hard factors. The `coordinates`
 * are interpreted as z, (y if applicable), (xiprime if applicable), rx, ry, etc.
 */
double gsl_monte_wrapper(double* coordinates, size_t ncoords, void* closure) {
    double real;
    // this does basically the same thing as cubature_wrapper but with a different signature
    cubature_wrapper(static_cast<unsigned int>(ncoords),  coordinates,  closure, 1, &real);
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
    if (*p_abserr != 0) {
        do {
            gsl_monte_vegas_integrate(&f, min, max, dim, incremental_iterations, rng, s, p_result, p_abserr);
            checkfinite(*p_result);
            checkfinite(*p_abserr);
            if (callback) {
                (*callback)(p_result, p_abserr, s);
            }
        } while (*p_abserr > 0 && fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.2);
    }
    gsl_monte_vegas_free(s);
    s = NULL;
}

/**
 * A wrapper function that calls the integration routine for quasi Monte Carlo integration.
 *
 * @param func the integrand
 * @param dim the number of dimensions it's being integrated over
 * @param closure something to be passed to the integrand as its last argument
 * @param min the lower bounds of the integration region
 * @param max the upper bounds of the integration region
 * @param[out] p_result the result
 * @param[out] p_abserr the error bound
 * @param iterations the maximum number of function evaluations
 * @param abserr the absolute error at which to stop
 * @param relerr the relative error at which to stop
 * @param rng the random number generator
 * @param callback a callback to call when the integration is done
 */
void quasi_integrate(double (*func)(double*, size_t, void*), size_t dim, void* closure, double* min, double* max, double* p_result, double* p_abserr,
               size_t iterations, double relerr, double abserr, gsl_qrng* qrng, void (*callback)(double*, double*, quasi_monte_state*)) {
    gsl_monte_function f;
    f.f = func;
    f.dim = dim;
    f.params = closure;

    quasi_monte_state* s = quasi_monte_alloc(dim);
    quasi_monte_integrate(&f, min, max, dim, iterations, relerr, abserr, qrng, s, p_result, p_abserr);
    checkfinite(*p_result);
    checkfinite(*p_abserr);
    if (callback) {
        (*callback)(p_result, p_abserr, s);
    }
    quasi_monte_free(s);
    s = NULL;
}

void cubature_integrate(integrand func, size_t dim, void* closure, double* min, double* max, double* p_result, double* p_abserr,
                        size_t iterations, double relerr, double abserr, void (*callback)(double*, double*)) {
    adapt_integrate(1, func, closure, static_cast<unsigned int>(dim), min, max, static_cast<unsigned int>(iterations), abserr, relerr, p_result, p_abserr);
    checkfinite(*p_result);
    checkfinite(*p_abserr);
    if (callback) {
        (*callback)(p_result, p_abserr);
    }
}

void Integrator::integrate_impl(double* result, double* error) {
    // it should already have been checked that there is at least one term of the appropriate type
    // and the type should be set appropriately
    size_t dimensions = current_integration_region->dimensions;
    double min[10];
    double max[10];
    assert(sizeof(min) / sizeof(min[0]) >= dimensions);
    assert(sizeof(max) / sizeof(max[0]) >= dimensions);
    current_integration_region->fill_min(ictx.ctx, min);
    current_integration_region->fill_max(ictx.ctx, max);
    switch (dimensions) {
        case 1:
        case 2:
            cubature_integrate(cubature_wrapper, dimensions, this, min, max, result, error, ictx.ctx.cubature_iterations, ictx.ctx.relerr, ictx.ctx.abserr, cubature_callback);
            break;
        default:
            if (ictx.ctx.strategy == MC_QUASI) {
                gsl_qrng* qrng = gsl_qrng_alloc(ictx.ctx.quasirandom_generator_type, static_cast<unsigned int>(dimensions));
                quasi_integrate(gsl_monte_wrapper, dimensions, this, min, max, result, error, ictx.ctx.quasi_iterations, ictx.ctx.relerr, ictx.ctx.abserr, qrng, quasi_callback);
                gsl_qrng_free(qrng);
                qrng = NULL;
            }
            else {
                gsl_rng* rng = gsl_rng_alloc(ictx.ctx.pseudorandom_generator_type);
                gsl_rng_set(rng, ictx.ctx.pseudorandom_generator_seed);
                switch (ictx.ctx.strategy) {
                    case MC_VEGAS:
                        vegas_integrate(gsl_monte_wrapper, dimensions, this, min, max, result, error, ictx.ctx.vegas_initial_iterations, ictx.ctx.vegas_incremental_iterations, rng, vegas_callback);
                        break;
                    case MC_MISER:
                        miser_integrate(gsl_monte_wrapper, dimensions, this, min, max, result, error, ictx.ctx.miser_iterations, rng, miser_callback);
                        break;
                    case MC_PLAIN:
                        throw "Unsupported integration method PLAIN";
                    default:
                        throw "Unknown integration method";
                }
                gsl_rng_free(rng);
                rng = NULL;
            }
    }
    if (callback) {
        callback(NULL, 0, 0);
    }
}

void Integrator::integrate(double* real, double* imag, double* error) {
    double result = 0.0;
    double abserr = 0.0;
    double tmp_result, tmp_error;

    for (HardFactorTypeMap::iterator it = terms.begin(); it != terms.end(); it++) {
        assert(it->second.size() > 0);
        set_current_integration_type(it->first);

        xi_preintegrated_term = false;
        integrate_impl(&tmp_result, &tmp_error);
        result += tmp_result;
        abserr += tmp_error;

        xi_preintegrated_term = true;
        integrate_impl(&tmp_result, &tmp_error);
        result += tmp_result;
        abserr += tmp_error;
    }

    *real = result;
    *imag = 0;
    *error = abserr;
}
