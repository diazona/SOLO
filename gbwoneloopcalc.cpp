#include <cassert>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_gamma.h>
#ifdef MONTECARLO
#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#endif
#include <iostream>
#include <cstdlib>
#include "cubature.h"
#include "mstwpdf.h"
#include "dss_pinlo.h"

using namespace std;

const int SUCCESS = 0;

/** Euler-Mascheroni constant. Value is copy-pasted from Wikipedia. */
const double EULER_GAMMA = 0.57721566490153286060651209008240243104215933593992;

void eprint(char* s) { cerr << s << endl; }

class Context {
public:
    double x0;
    double A;
    double c;
    double lambda;
    double mu2;
    double Nc;
    double Nf;
    double Sperp;
    double pT2;
    double sqs;
    double Y;
    double alphasbar_fixed;

    c_mstwpdf* pdf_object;
    DSSpiNLO* ff_object;
    
    Context(const char* pdf_filename, const char* ff_filename) :
     x0(0.000304), A(1), c(1), lambda(0.288), mu2(1), Nc(3), Nf(3), Sperp(1), pT2(0), sqs(0), Y(0), alphasbar_fixed(0.1 / (2 * M_PI)) {
         ff_object = new DSSpiNLO(ff_filename);
         pdf_object = new c_mstwpdf(pdf_filename);
    }
    
    Context(double x0, double A, double c, double lambda, double mu2, double Nc, double Nf, double Sperp, double pT2, double sqs, double Y, double alphasbar_fixed, const char* pdf_filename, const char* ff_filename) :
     x0(x0), A(A), c(c), lambda(lambda), mu2(mu2), Nc(Nc), Nf(Nf), Sperp(Sperp), pT2(pT2), sqs(sqs), Y(Y), alphasbar_fixed(alphasbar_fixed) {
         ff_object = new DSSpiNLO(ff_filename);
         pdf_object = new c_mstwpdf(pdf_filename);
    }
    
    ~Context() {
        delete pdf_object;
        delete ff_object;
    }
};



/** A collection of global parameters for the calculation. */
class HardFactor {
public:
    // These fields should be considered read-only once set by the constructor
    double x0;
    double A;
    double c;
    double lambda;
    double Q02x0lambda;
    double mu2;
    double Nc;
    double Nf;
    double Sperp;
    double pT2;
    double sqs;
    double Y;
    double tau;
    double alphasbar_fixed;

    c_mstwpdf* pdf_object;
    DSSpiNLO* ff_object;
    
    Context* ctx;
    
    HardFactor(Context* ctx) :
     x0(ctx->x0),
     A(ctx->A),
     c(ctx->c),
     lambda(ctx->lambda),
     mu2(ctx->mu2),
     Nc(ctx->Nc),
     Nf(ctx->Nf),
     Sperp(ctx->Sperp),
     pT2(ctx->pT2),
     sqs(ctx->sqs),
     Y(ctx->Y),
     alphasbar_fixed(ctx->alphasbar_fixed),
     pdf_object(ctx->pdf_object),
     ff_object(ctx->ff_object),
     ctx(ctx) {
        Q02x0lambda = c * pow(A, 1.0d/3.0d) * pow(x0, lambda);
        tau = sqrt(pT2)/sqs*exp(Y);
    }

    double quarkfactor(double z, double xi);
    
    double hardFactorIntegral();
    double hard_integrand_1D(double z);
    double hard_integrand_2D(double z, double y);
    
    virtual double alphasbar(double q2) {
        return alphasbar_fixed; // fixed coupling
    }
    
    virtual double Fs(double, double) = 0;
    virtual double Fn(double, double) = 0;
    virtual double Fd(double, double) = 0;
    
    // saturation scale: Qs2(z, xi)
    virtual double Qs2(double z, double xi) {
        double xg = sqrt(pT2) / (sqs * z) * exp(-Y);
        return Q02x0lambda / pow(xg, lambda); // Q_0^2 (x_0 / x)^lambda
    }

    // gluon distribution: Fxg(q2, z, xi)
    virtual double Fxg(double q2, double z, double xi) = 0;
};

class GBWHardFactor : public HardFactor {
public:
    GBWHardFactor(Context* ctx) :
        HardFactor(ctx) {}
    double Fxg(double q2, double z, double xi) {
        double Qs2 = this->Qs2(z, xi);
        return Sperp / (M_PI * Qs2) * exp(-q2 / Qs2);
    }
};

// Integrating the hard factor

// the integration region is a square in the z-y plane
double HardFactor::hard_integrand_2D(double z, double y) {
    double result = 0.0d;
    double xi = (y * (z - tau) - tau * (z - 1)) / (z * (1 - tau));
    result += (Fs(z, xi) - Fs(z, 1)) / (1 - xi);
    result += Fn(z, xi);
    result *= (1 - tau / z) / (1 - tau); // Jacobian from y to xi
    return result;
}

double HardFactor::hard_integrand_1D(double z){
    double result = 0.0d;
    result += Fs(z, 1)*log(1 - tau/z);
    result += Fd(z, 1);
    return result;
}

#ifdef MONTECARLO
double gsl_monte_wrapper_1D(double* coordinates, size_t ncoords, void* closure) {
    assert(ncoords == 1);
    double value = ((HardFactor*)closure)->hard_integrand_1D(coordinates[0]);
    return value;
}

double gsl_monte_wrapper_2D(double* coordinates, size_t ncoords, void* closure) {
    assert(ncoords == 2);
    double value = ((HardFactor*)closure)->hard_integrand_2D(coordinates[0], coordinates[1]);
    return value;
}

void vegas_integrate(double (*func)(double*, size_t, void*), size_t dim, void* closure, double* min, double* max, double* p_result, double* p_abserr,
               void (*callback)(double*, double*, gsl_monte_vegas_state*)) {
    double old_result = NAN;
    gsl_monte_function f;
    f.f = func;
    f.dim = dim;
    f.params = closure;

    gsl_rng_env_setup();
    gsl_rng* rng = gsl_rng_alloc(gsl_rng_default);
    gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(dim);
    gsl_monte_vegas_integrate(&f, min, max, dim, 10000, rng, s, p_result, p_abserr);
    if (callback) {
        (*callback)(p_result, p_abserr, s);
    }
    do {
        gsl_monte_vegas_integrate(&f, min, max, dim, 100000, rng, s, p_result, p_abserr);
        if (callback) {
            (*callback)(p_result, p_abserr, s);
        }
        old_result = *p_result;
    } while (old_result != *p_result && fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.2);
    gsl_monte_vegas_free(s);
    s = NULL;
    gsl_rng_free(rng);
    rng = NULL;
}

void eprint_callback(double* p_result, double* p_abserr, gsl_monte_vegas_state* s) {
    cerr << "VEGAS output: " << *p_result << " err: " << *p_abserr << " chisq:" << gsl_monte_vegas_chisq(s) << endl;
}
#else
void cubature_wrapper_1D(unsigned int ncoords, const double* coordinates, void* closure, unsigned int nvalues, double* values) {
    assert(ncoords == 1);
    assert(nvalues == 1);
    values[0] = ((HardFactor*)closure)->hard_integrand_1D(coordinates[0]);
}

void cubature_wrapper_2D(unsigned int ncoords, const double* coordinates, void* closure, unsigned int nvalues, double* values) {
    assert(ncoords == 2);
    assert(nvalues == 1);
    values[0] = ((HardFactor*)closure)->hard_integrand_2D(coordinates[0], coordinates[1]);
}
#endif

/** Computes the integral of a hard factor over z and xi. */
double HardFactor::hardFactorIntegral() {
    double result = 0.0d;
    double tmp_result = 0.0d;
    double tmp_error = 0.0d;
    double xmin1D[] = {tau};
    double xmax1D[] = {1.0d};
    double xmin2D[] = {tau, tau};
    double xmax2D[] = {1.0d, 1.0d};
    int status = 0;
#ifdef MONTECARLO
    // the 2D integration
    vegas_integrate(gsl_monte_wrapper_2D, 2, this, xmin2D, xmax2D, &tmp_result, &tmp_error, eprint_callback);
    result += tmp_result;
    // the 1D integration
    vegas_integrate(gsl_monte_wrapper_1D, 1, this, xmin1D, xmax1D, &tmp_result, &tmp_error, eprint_callback);
    result += tmp_result;
#else
    // the 2D integration
    status = adapt_integrate(1, cubature_wrapper_2D, this, 2, xmin2D, xmax2D, 100000, 0, 1e-3, &tmp_result, &tmp_error);
    if (status != SUCCESS) {
        cerr << "Error in 2D integration (probably memory)" << endl;
        exit(1);
    }
    result += tmp_result;
    // the 1D integration
    status = adapt_integrate(1, cubature_wrapper_1D, this, 1, xmin1D, xmax1D, 100000, 0, 1e-3, &tmp_result, &tmp_error);
    if (status != SUCCESS) {
        cerr << "Error in 1D integration (probably memory)" << endl;
        exit(1);
    }
    result += tmp_result;
#endif
    return result;
}

double HardFactor::quarkfactor(double z, double xi) {
    double x = tau / (z * xi);
    double result = 0.0d;
    pdf_object->update(x, sqrt(mu2));
    ff_object->update(z, mu2);
    
    // Proton contributions:
    // up quark
    result += (pdf_object->cont.upv + pdf_object->cont.usea) * ff_object->fragmentation(DSSpiNLO::up, DSSpiNLO::pi_minus);
    // antiup quark
    result += pdf_object->cont.usea * ff_object->fragmentation(DSSpiNLO::up_bar, DSSpiNLO::pi_minus);
    // down quark
    result += (pdf_object->cont.dnv + pdf_object->cont.dsea) * ff_object->fragmentation(DSSpiNLO::down, DSSpiNLO::pi_minus);
    // antidown quark
    result += pdf_object->cont.dsea * ff_object->fragmentation(DSSpiNLO::down_bar, DSSpiNLO::pi_minus);
    // strange quark
    result += pdf_object->cont.str * ff_object->fragmentation(DSSpiNLO::strange, DSSpiNLO::pi_minus);
    // antistrange quark
    result += pdf_object->cont.sbar * ff_object->fragmentation(DSSpiNLO::strange_bar, DSSpiNLO::pi_minus);
    
    // Neutron contributions (for deuteron collisions), assuming isospin symmetry:
    // up quark
    result += (pdf_object->cont.dnv + pdf_object->cont.dsea) * ff_object->fragmentation(DSSpiNLO::up, DSSpiNLO::pi_minus);
    // antiup quark
    result += pdf_object->cont.dsea * ff_object->fragmentation(DSSpiNLO::up_bar, DSSpiNLO::pi_minus);
    // down quark
    result += (pdf_object->cont.upv + pdf_object->cont.usea) * ff_object->fragmentation(DSSpiNLO::down, DSSpiNLO::pi_minus);
    // antidown quark
    result += pdf_object->cont.usea * ff_object->fragmentation(DSSpiNLO::down_bar, DSSpiNLO::pi_minus);
    // strange quark
    result += pdf_object->cont.str * ff_object->fragmentation(DSSpiNLO::strange, DSSpiNLO::pi_minus);
    // antistrange quark
    result += pdf_object->cont.sbar * ff_object->fragmentation(DSSpiNLO::strange_bar, DSSpiNLO::pi_minus);
    
    return result;
}

// Calculation of H02qq

class H02qq : public GBWHardFactor {
public:
    H02qq(Context* ctx) :
        GBWHardFactor(ctx) {}
    double Fs(double, double) { return 0.0d; };
    double Fn(double, double) { return 0.0d; };
    double Fd(double, double);
};

double H02qq::Fd(double z, double xi) {
    double kT2 = pT2/(z*z);
    double Qs2 = this->Qs2(z, xi);
    return Sperp / (M_PI * Qs2) * quarkfactor(z, xi) / (z * z) * exp(-kT2 / Qs2);
}

// Calculation of H12qq

class H12qq : public GBWHardFactor {
public:
    H12qq(Context* ctx) :
        GBWHardFactor(ctx) {}
    double Fs(double, double);
    double Fn(double, double) { return 0.0d; };
    double Fd(double, double);
};

/**
 * Returns the product of the exponent and the multivariate Laguerre polynomial.
 *  exp(x) L'(-1, -x)
 * or equivalently
 *  -gamma_E - Gamma(0, -x) - ln(-x)
 * or equivalently (and best for implementation if x > 0)
 *  -gamma_E + Ei(x) - ln(x)
 */
inline double exp_multi_laguerre(double x) {
    if (x < 0) {
        return -EULER_GAMMA - gsl_sf_gamma_inc(0, -x) - log(-x);
    }
    else {
        return -EULER_GAMMA + gsl_sf_expint_Ei(x) - log(x);
    }
}

inline double h12qq_singularterm(double xi, double kT2, double Qs2, HardFactor* hf) {
    double r2 = kT2 / Qs2;
    return 0.5 * hf->Nc                     // Nc/2
     * hf->Sperp / (M_PI * Qs2) * exp(-r2)  // Fg(kT)
     * (log(Qs2 / hf->mu2) - EULER_GAMMA    // ln(Qs2/mu2 e^gamma_E)
            + exp_multi_laguerre(r2));      // exp(kT2/Qs2)L'(-1, -kT2/Qs2)
}

double H12qq::Fs(double z, double xi) {
    double kT2 = pT2/(z*z);
    double Qs2 = this->Qs2(z, xi);
    double xi2 = xi * xi;
    return alphasbar(kT2) * quarkfactor(z, xi) / (z*z) * (1 + xi2) * (
        h12qq_singularterm(xi, kT2, Qs2, this)
      + h12qq_singularterm(xi, kT2 / xi2, Qs2, this) / xi2
    );
}

double H12qq::Fd(double z, double xi) {
    double kT2 = pT2/(z*z);
    double Qs2 = this->Qs2(z, xi);
    double r2 = kT2 / Qs2;
    double xi2 = xi * xi;
    return alphasbar(kT2) * quarkfactor(z, xi) / (z*z) * 1.5 * (
        h12qq_singularterm(xi, kT2, Qs2, this)
      + h12qq_singularterm(xi, kT2 / xi2, Qs2, this) / xi2
      - Nc
        * Sperp / (M_PI * Qs2) * exp(-r2) // Fg(kT)
        * (log(Qs2 / kT2) - EULER_GAMMA   // ln(Qs2/kT2 e^gamma_E)
            + exp_multi_laguerre(r2))     // exp(kT2/Qs2)L'(-1, -kT2/Qs2)
    );
}

// Calculation of H14

class H14qq : public GBWHardFactor {
public:
    H14qq(Context* ctx) :
        GBWHardFactor(ctx) {}
    double Fs(double, double);
    double Fn(double, double) { return 0.0d; };
    double Fd(double, double);
};

double H14qq::Fs(double z, double xi) {
    double kT2 = pT2/(z*z);
    double Qs2 = this->Qs2(z, xi);
    return -alphasbar(kT2) * quarkfactor(z, xi) / (z*z) * Nc * Sperp / M_PI * (1 + xi*xi) / kT2
        * (1 - exp(-kT2 / Qs2)) * (1 - exp(-kT2 / Qs2 / (xi*xi)));
}

void cubature_wrapper_internal_integrand(unsigned int ncoords, const double* coordinates, void* closure, unsigned int nvalues, double* values) {
    assert(ncoords == 1);
    assert(nvalues == 1);
    double xip = coordinates[0];
    double xip2 = xip * xip;
    double x2 = *((double*)closure);
    values[0] = ((1 + xip2) * exp_multi_laguerre(-xip2 * x2) - 2 * exp_multi_laguerre(-x2)) / (1 - xip);
}

double h14qq_internalIntegral(double x2) {
    double result = 0.0d;
    double error = 0.0d;
    double xmin1D[] = {0.0d};
    double xmax1D[] = {1.0d};
    int status = adapt_integrate(1, cubature_wrapper_internal_integrand, &x2, 1, xmin1D, xmax1D, 10000, 0, 1e-4, &result, &error);
    if (status != SUCCESS) {
        cerr << "Error in 1D integration (probably memory)" << endl;
        exit(1);
    }
    return result;
}

double H14qq::Fd(double z, double xi) {
    double kT2 = pT2/(z*z);
    double Qs2 = this->Qs2(z, xi);
    double r2 = kT2 / Qs2;
    return alphasbar(kT2) * quarkfactor(z, xi) / (z*z) * Nc 
        * Sperp / M_PI / Qs2 * exp(-r2)       // Fg(kT)
        * (-1.5 * (log(r2) + EULER_GAMMA) 
            + h14qq_internalIntegral(r2));
}

double calculateLOterm(Context* ctx) {
    H02qq* h02qq = new H02qq(ctx);
    double result = h02qq->hardFactorIntegral();
    delete h02qq;
    return result;
}

double calculateNLOterm(Context* ctx) {
    H12qq* h12qq = new H12qq(ctx);
    H14qq* h14qq = new H14qq(ctx);
    double result12 = h12qq->hardFactorIntegral();
    double result14 = h14qq->hardFactorIntegral();
    delete h12qq;
    delete h14qq;
    return result12 + result14;
}

void fillYieldArray(double sqs, double Y, int pTlen, double* pT, double* yield) {
    int i;
    Context gctx("mstw2008nlo.00.dat", "PINLO.DAT");
    gctx.Sperp = 1.0d;
    gctx.mu2 = 10;
    gctx.Nc = 3;
    gctx.Nf = 3;
    gctx.sqs = sqs;
    gctx.Y = Y;
    gctx.A = 197;
    gctx.c = 0.56;
    gctx.x0 = 0.000304;
    gctx.lambda = 0.288;
    gctx.alphasbar_fixed = 0.2 / (2 * M_PI);
    for (i = 0; i < pTlen; i++) {
        gctx.pT2 = pT[i]*pT[i];
        yield[2*i] = calculateLOterm(&gctx);
        yield[2*i+1] = calculateNLOterm(&gctx);
    }
}

int main(int argc, char** argv) {
    double pT[] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4};
    size_t len = sizeof(pT)/sizeof(double);
    double yield[2*len];
    
    fillYieldArray(200, 3.2, len, pT, yield);
    cout << "pT\tLOqq\tNLOqq\tLOqq+NLOqq" << endl;
    for (size_t i = 0; i < len; i++) {
        cout << pT[i] << "\t" << yield[2*i] << "\t" << yield[2*i+1] << "\t" << yield[2*i] + yield[2*i+1] << endl;
    }
    return 0;
}

