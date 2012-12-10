#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include "context.h"
#include "integrationcontext.h"
#include "hardfactor.h"

static enum {MC_PLAIN, MC_MISER, MC_VEGAS} integration_strategy = MC_VEGAS;

class Integrator {
private:
    IntegrationContext* ictx;
    HardFactor** dipole_terms;
    size_t n_dipole_terms;
    HardFactor** quadrupole_terms;
    size_t n_quadrupole_terms;
    HardFactor** momentum1_terms;
    size_t n_momentum1_terms;
    HardFactor** momentum2_terms;
    size_t n_momentum2_terms;
    HardFactor** momentum3_terms;
    size_t n_momentum3_terms;
    term_type current_term_type;
    void (*callback)(IntegrationContext*, double, double);
public:
    Integrator(Context* ctx, size_t hflen, HardFactor** hflist);
    ~Integrator();
    void update_position(double z, double y, double rx, double ry);
    void update_position(double z, double y, double sx, double sy, double tx, double ty);
    void update_momentum(double z, double y, double q1x, double q1y);
    void update_momentum(double z, double y, double q1x, double q1y, double q2x, double q2y);
    void update_momentum(double z, double y, double q1x, double q1y, double q2x, double q2y, double q3x, double q3y);
    void evaluate_1D_integrand(double* real, double* imag);
    void evaluate_2D_integrand(double* real, double* imag);
    void integrate(double* real, double* imag, double* error);
    void set_current_term_type(term_type new_term_type) {
        current_term_type = new_term_type;
    }
    void set_callback(void (*callback)(IntegrationContext*, double, double)) {
        this->callback = callback;
    }
private:
    void integrate_impl(double (*monte_wrapper)(double*, size_t, void*), size_t dimensions, double* min, double* max, double* result, double* error);
};

#endif // _INTEGRATOR_H_
