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

#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include <map>
#include <vector>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include "../configuration/context.h"
#include "integrationcontext.h"
#include "integrationtype.h"
#include "../hardfactors/hardfactor.h"
#include "quasimontecarlo.h"
#include "../typedefs.h"

/**
 * A class to interface with the GSL Monte Carlo integration routines.
 * 
 * Integrator objects basically have two purposes: most importantly,
 * the GSL integration functions accept a void pointer that is passed on to
 * the function being evaluated, which can be used to pass an object that
 * contains all the extra information needed to evaluate the function.
 * For this program, the object passed is an instance of Integrator. It
 * contains the map of hard factors and the IntegrationContext.
 * 
 * Integrator also contains the code to prepare for and actually invoke
 * the integration functions. It sorts the hard factors by integration type,
 * sets up the functions to be evaluated, creates the random number generator,
 * and so on.
 * 
 * For functions that have 1D and 2D versions: 1D refers to integrating
 * over z only, whereas 2D refers to integrating over z and y (equivalently,
 * z and xi). The integration that is actually performed will have many more
 * than 1 or 2 dimensions: the "1D" integral is (1+N)-dimensional and the
 * "2D" integral is (2+N)-dimensional, where N is the number of dimensions
 * specified by the current IntegrationType.
 */
class Integrator {
private:
    /**
     * The type of term being integrated.
     * 
     * Types are defined in integrationtype.h.
     */
    IntegrationType const* current_type;
    /**
     * The IntegrationContext object associated with this Integrator.
     */
    IntegrationContext ictx;
    /**
     * The hard factors that this Integrator should be integrating, sorted
     * out by their type.
     */
    HardFactorTypeMap terms;
    /** A callback function to call each time the function is evaluated */
    void (*callback)(const IntegrationContext*, double, double);
    /** A callback function to call each time a cubature integration finishes */
    void (*cubature_callback)(double*, double*);
    /** A callback function to call each time a MISER integration finishes */
    void (*miser_callback)(double*, double*, gsl_monte_miser_state*);
    /** A callback function to call each time a step of the VEGAS integration finishes */
    void (*vegas_callback)(double*, double*, gsl_monte_vegas_state*);
    /** A callback function to call each time a quasi Monte Carlo integration finishes */
    void (*quasi_callback)(double*, double*, quasi_monte_state*);

    size_t core_dimensions;
    
    double xg_min, xg_max;
public:
    Integrator(const Context* ctx, const ThreadLocalContext* tlctx, const HardFactorList& hflist, const double xg_min, const double xg_max);
    ~Integrator();
    /**
     * Evaluates the 1D or 2D integrand, using the values of variables currently stored
     * in the IntegrationContext, and stores the real and imaginary parts of the
     * result in the given variables.
     */
    void evaluate_integrand(double* real, double* imag);
    /**
     * Performs the 1D and 2D integrals, and stores the result and error bound
     * in the variables `real` and `error`. The variable `imag` is set to zero
     * in this implementation because it's theoretically supposed to be zero,
     * but it's in the function declaration in case we wanted to calculate
     * the imaginary part for some reason.
     * 
     * What is actually integrated by this function is the set of all hard factors
     * whose IntegrationType matches the current one. The current IntegrationType
     * is set by the method set_current_integration_type(). So the overall workflow
     * is to call set_current_integration_type() and then integrate() for each type
     * in turn.
     */
    void integrate(double* real, double* imag, double* error);
    /**
     * Sets the current integration type.
     */
    void set_current_integration_type(IntegrationType const* new_type) {
        current_type = new_type;
    }
    /**
     * Sets the callback to be invoked each time the integrand is evaluated.
     */
    void set_callback(void (*callback)(const IntegrationContext*, double, double)) {
        this->callback = callback;
    }
    /**
     * Sets the callback to be invoked each time a MISER integration finishes.
     */
    void set_cubature_callback(void (*cubature_callback)(double*, double*)) {
        this->cubature_callback = cubature_callback;
    }
    /**
     * Sets the callback to be invoked each time a MISER integration finishes.
     */
    void set_miser_callback(void (*miser_callback)(double*, double*, gsl_monte_miser_state*)) {
        this->miser_callback = miser_callback;
    }
    /**
     * Sets the callback to be invoked each time a step of a VEGAS integration finishes.
     */
    void set_vegas_callback(void (*vegas_callback)(double*, double*, gsl_monte_vegas_state*)) {
        this->vegas_callback = vegas_callback;
    }
    /**
     * Sets the callback to be invoked each time a quasi Monte Carlo integration finishes.
     */
    void set_quasi_callback(void (*quasi_callback)(double*, double*, quasi_monte_state*)) {
        this->quasi_callback = quasi_callback;
    }
private:
    /**
     * Most of the code for integrating "1D" and "2D" integrals is the same, so it's
     * abstracted into this method.
     */
    void integrate_impl(const size_t core_dimensions, double* result, double* error);
    
    // the non-member functions that actually implement the integration
    friend void cubature_wrapper(unsigned int ncoords, const double* coordinates, void* closure, unsigned int nresults, double* results);
    friend double gsl_monte_wrapper(double* coordinates, size_t ncoords, void* closure);
};

#endif // _INTEGRATOR_H_
