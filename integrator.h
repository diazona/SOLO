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
    term_type current_term_type;
    void (*callback)(IntegrationContext*, double, double);
public:
    Integrator(Context* ctx, size_t hflen, HardFactor** hflist);
    ~Integrator();
    void update(double z, double y, double rx, double ry) {
        // I need (z = z, y = y, xx = sx + bx, xy = sy + by, yx = tx + bx, yy = ty + by)
        // we assume no impact parameter dependence, so I can arbitrarily set bx, by, tx, ty to 0
        assert(current_term_type == dipole);
        ictx->update(z, y, rx, ry, 0, 0, 0, 0);
    }
    void update(double z, double y, double sx, double sy, double tx, double ty) {
        // I need (z = z, y = y, xx = sx + bx, xy = sy + by, yx = tx + bx, yy = ty + by)
        // we assume no impact parameter dependence, so I can arbitrarily set bx, by to 0
        assert(current_term_type == quadrupole);
        ictx->update(z, y, sx, sy, tx, ty, 0, 0);
    }
    void evaluate_1D_integrand(double* real, double* imag);
    void evaluate_2D_integrand(double* real, double* imag);
    void integrate(double* real, double* imag, double* error);
    void set_current_term_type(term_type new_term_type) {
        current_term_type = new_term_type;
    }
    void set_callback(void (*callback)(IntegrationContext*, double, double)) {
        this->callback = callback;
    }
};

#endif // _INTEGRATOR_H_
