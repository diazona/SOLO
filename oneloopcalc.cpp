#include <cassert>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include "cubature.h"
#include "mstwpdf.h"
#include "dss_pinlo.h"

using namespace std;

const int SUCCESS = 0;

static enum {MC_PLAIN, MC_MISER, MC_VEGAS} integration_strategy = MC_VEGAS;
typedef enum {NONE=0, dipole=2, quadrupole=4} term_type;
static bool trace = false;

/** Euler-Mascheroni constant. Value is copy-pasted from Wikipedia. */
const double EULER_GAMMA = 0.57721566490153286060651209008240243104215933593992;

void eprint(char* s) { cerr << s << endl; }

class IntegrationContext;
class GluonDistribution;
class HardFactor;

class Context {
public:
    double x0;
    double A;
    double c;
    double lambda;
    double mu2;
    double Nc;
    double Nf;
    double CF;
    double TR;
    double Sperp;
    double pT2;
    double sqs;
    double Y;
    double Q02x0lambda;
    double tau;
    double alphasbar_fixed;
    
    GluonDistribution* gdist;

    // NOTE: these contain state that should be associated with IntegrationContext.
    // So don't use one Context with more than one IntegrationContext at once.
    c_mstwpdf* pdf_object;
    DSSpiNLO* ff_object;

    Context(double x0, double A, double c, double lambda, double mu2, double Nc, double Nf, double CF, double TR, double Sperp, double pT2, double sqs, double Y, double alphasbar_fixed, GluonDistribution* gdist, const char* pdf_filename, const char* ff_filename) :
     x0(x0), A(A), c(c), lambda(lambda), mu2(mu2), Nc(Nc), Nf(Nf), CF(CF), TR(TR), Sperp(Sperp), pT2(pT2), sqs(sqs), Y(Y), alphasbar_fixed(alphasbar_fixed), gdist(gdist) {
         recalculate();
         pdf_object = new c_mstwpdf(pdf_filename);
         ff_object = new DSSpiNLO(ff_filename);
    }

    ~Context() {
        if (pdf_object) {
            delete pdf_object;
        }
        if (ff_object) {
            delete ff_object;
        }
    }

    void recalculate() {
         Q02x0lambda = c * pow(A, 1.0d/3.0d) * pow(x0, lambda);
         tau = sqrt(pT2)/sqs*exp(Y);
    }
    
    virtual double alphasbar(double kT2) {
        return alphasbar_fixed;
    }
};

class GluonDistribution {
public:
    virtual double S2(double r2, IntegrationContext* ictx) = 0;
    virtual double S4(double r2, double s2, double t2, IntegrationContext* ictx) = 0;
};

class GBWGluonDistribution: public GluonDistribution {
public:
    double S2(double r2, IntegrationContext* ictx);
    double S4(double r2, double s2, double t2, IntegrationContext* ictx);
};

class MVGluonDistribution: public GluonDistribution {
public:
    MVGluonDistribution(double LambdaMV) : LambdaMV(LambdaMV) {};
    double S2(double r2, IntegrationContext* ictx);
    double S4(double r2, double s2, double t2, IntegrationContext* ictx);
private:
    double LambdaMV;
};

// Not even anywhere close to thread-safe!
class IntegrationContext {
public:
    Context* ctx;
    // updated
    double z;
    double xi;
    double xx, xy;
    double yx, yy;
    double bx, by;
    // calculated
    double z2, xi2;
    double kT2, kT;
    double xp, xg;
    double r2, s2, t2;
    double Qs2;
    double alphasbar;
    double qqfactor;
    double ggfactor;
    double gqfactor;
    double S2r, S4rst;
    
    IntegrationContext(Context* ctx) :
      ctx(ctx),
      z(0), xi(0),
      xx(0), xy(0),
      yx(0), yy(0),
      bx(0), by(0),
      z2(0), xi2(0),
      kT2(0), kT(0),
      xp(0), xg(0),
      r2(0), s2(0), t2(0),
      Qs2(0),
      alphasbar(0),
      qqfactor(0),
      ggfactor(0),
      gqfactor(0),
      S2r(0), S4rst(0) {
    };
    void update(double z, double y, double xx, double xy, double yx, double yy, double bx, double by);
};

void IntegrationContext::update(double z, double y, double xx, double xy, double yx, double yy, double bx, double by) {
    double qqfactor = 0.0d, ggfactor = 0.0d, gqfactor = 0.0d;
    c_mstwpdf* pdf_object = ctx->pdf_object;
    DSSpiNLO* ff_object = ctx->ff_object;
    assert(z <= 1);
    assert(z >= ctx->tau);
    assert(y <= 1);
    assert(y >= ctx->tau);
    assert(xx == xx);
    assert(xy == xy);
    assert(yx == yx);
    assert(yy == yy);
    assert(bx == bx);
    assert(by == by);
    this->z = z;
    this->z2 = z*z;
    if (y == 1.0d) {
        this->xi = 1.0d; // avoid floating-point roundoff error
    }
    else {
        this->xi = (y * (z - ctx->tau) - ctx->tau * (z - 1)) / (z * (1 - ctx->tau));
    }
    this->xi2 = xi*xi;
    this->xx = xx;
    this->xy = xy;
    this->yx = yx;
    this->yy = yy;
    this->bx = bx;
    this->by = by;
    this->r2 = (xx - yx) * (xx- yx) + (xy - yy) * (xy - yy);
    this->s2 = (xx - bx) * (xx- bx) + (xy - by) * (xy - by);
    this->t2 = (yx - bx) * (yx- bx) + (yy - by) * (yy - by);
    this->xp = ctx->tau / (z * xi);
    this->kT2 = ctx->pT2 / this->z2;
    this->kT = sqrt(this->kT2);
    this->xg = kT / ctx->sqs * exp(-ctx->Y);
    this->Qs2 = ctx->Q02x0lambda / pow(this->xg, ctx->lambda); // Q_0^2 (x_0 / x)^lambda
    this->alphasbar = ctx->alphasbar(this->kT2);

    // Calculate the new gluon distribution values
    // this has to be done after kinematics are updated
    this->S2r = ctx->gdist->S2(r2, this);
    this->S4rst = ctx->gdist->S4(r2, s2, t2, this);

    // Calculate the new quark/gluon factors
    pdf_object->update(xp, sqrt(ctx->mu2));
    ff_object->update(z, ctx->mu2);
    
    // Proton contributions:
    qqfactor += (pdf_object->cont.upv + pdf_object->cont.usea) * ff_object->fragmentation(DSSpiNLO::up, DSSpiNLO::pi_minus);
    qqfactor += pdf_object->cont.usea * ff_object->fragmentation(DSSpiNLO::up_bar, DSSpiNLO::pi_minus);
    qqfactor += (pdf_object->cont.dnv + pdf_object->cont.dsea) * ff_object->fragmentation(DSSpiNLO::down, DSSpiNLO::pi_minus);
    qqfactor += pdf_object->cont.dsea * ff_object->fragmentation(DSSpiNLO::down_bar, DSSpiNLO::pi_minus);
    qqfactor += pdf_object->cont.str * ff_object->fragmentation(DSSpiNLO::strange, DSSpiNLO::pi_minus);
    qqfactor += pdf_object->cont.sbar * ff_object->fragmentation(DSSpiNLO::strange_bar, DSSpiNLO::pi_minus);
    
    // Neutron contributions (for deuteron collisions), assuming isospin symmetry:
    qqfactor += (pdf_object->cont.dnv + pdf_object->cont.dsea) * ff_object->fragmentation(DSSpiNLO::up, DSSpiNLO::pi_minus);
    qqfactor += pdf_object->cont.dsea * ff_object->fragmentation(DSSpiNLO::up_bar, DSSpiNLO::pi_minus);
    qqfactor += (pdf_object->cont.upv + pdf_object->cont.usea) * ff_object->fragmentation(DSSpiNLO::down, DSSpiNLO::pi_minus);
    qqfactor += pdf_object->cont.usea * ff_object->fragmentation(DSSpiNLO::down_bar, DSSpiNLO::pi_minus);
    qqfactor += pdf_object->cont.str * ff_object->fragmentation(DSSpiNLO::strange, DSSpiNLO::pi_minus);
    qqfactor += pdf_object->cont.sbar * ff_object->fragmentation(DSSpiNLO::strange_bar, DSSpiNLO::pi_minus);
    
    this->qqfactor = qqfactor;

    
    // Proton contribution:
    ggfactor += pdf_object->cont.glu * ff_object->fragmentation(DSSpiNLO::gluon, DSSpiNLO::pi_minus);
    
    // Neutron contribution (for deuteron collisions), assuming isospin symmetry:
    ggfactor *= 2;
    
    this->ggfactor = ggfactor;

    
    // Proton contributions:
    gqfactor += (  pdf_object->cont.upv + 2*pdf_object->cont.usea
                 + pdf_object->cont.dnv + 2*pdf_object->cont.dsea
                 + pdf_object->cont.str
                 + pdf_object->cont.sbar
                ) * ff_object->fragmentation(DSSpiNLO::gluon, DSSpiNLO::pi_minus);
    
    // Neutron contributions (for deuteron collisions), assuming isospin symmetry:
    gqfactor *= 2;
    
    this->gqfactor = gqfactor;
}

double GBWGluonDistribution::S2(double r2, IntegrationContext* ictx) {
    return exp(-0.25 * r2 * ictx->Qs2);
}
double GBWGluonDistribution::S4(double r2, double s2, double t2, IntegrationContext* ictx) {
    return exp(-0.25 * ictx->Qs2 * (s2 + t2));
}

double MVGluonDistribution::S2(double r2, IntegrationContext* ictx) {
    return pow(M_E + 1.0 / (sqrt(r2) * LambdaMV), -0.25 * r2 * ictx->Qs2);
}
double MVGluonDistribution::S4(double r2, double s2, double t2, IntegrationContext* ictx) {
    return S2(s2, ictx) * S2(t2, ictx);
}

class HardFactor {
public:
    virtual term_type get_type() = 0;
    virtual double real_singular_contribution(IntegrationContext* ictx) { return 0; }
    virtual double imag_singular_contribution(IntegrationContext* ictx) { return 0; }
    virtual double real_normal_contribution(IntegrationContext* ictx) { return 0; }
    virtual double imag_normal_contribution(IntegrationContext* ictx) { return 0; }
    virtual double real_delta_contribution(IntegrationContext* ictx) { return 0; }
    virtual double imag_delta_contribution(IntegrationContext* ictx) { return 0; }
};

class H02qq : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    double real_delta_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI) * ictx->qqfactor / ictx->z2 * ictx->S2r *
         // real part of the hard factor
         cos(ictx->kT * (ictx->xx - ictx->yx)); // take angle of k_perp to be 0
    }
    double imag_delta_contribution(IntegrationContext* ictx) {
        return -1/(4*M_PI*M_PI) * ictx->qqfactor / ictx->z2 * ictx->S2r *
         // imaginary part of the hard factor
         sin(ictx->kT * (ictx->xx - ictx->yx)); 
    }
};

class H12qq : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    double real_singular_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->qqfactor / ictx->z2 * ictx->S2r *
         // real part of the singular contribution
         ictx->ctx->CF * (1 + ictx->xi2) * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
         * (cos(ictx->kT * (ictx->xx - ictx->yx)) + cos(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi) / ictx->xi2);
    }
    double imag_singular_contribution(IntegrationContext* ictx) {
        return -1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->qqfactor / ictx->z2 * ictx->S2r *
         // imaginary part of the singular contribution
         ictx->ctx->CF * (1 + ictx->xi2) * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
         * (sin(ictx->kT * (ictx->xx - ictx->yx)) + sin(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi) / ictx->xi2);
    }
    double real_delta_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->qqfactor / ictx->z2 * ictx->S2r * (
         // real part of the delta contribution
         1.5 * ictx->ctx->CF * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
         * (cos(ictx->kT * (ictx->xx - ictx->yx)) + cos(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi) / ictx->xi2)
         - 3.0 * ictx->ctx->CF * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->kT2))
         * cos(ictx->kT * (ictx->xx - ictx->yx))
        );
    }
    double imag_delta_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->qqfactor / ictx->z2 * ictx->S2r * (
         // imaginary part of the delta contribution
         1.5 * ictx->ctx->CF * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
         * (sin(ictx->kT * (ictx->xx - ictx->yx)) + sin(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi) / ictx->xi2)
         - 3.0 * ictx->ctx->CF * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->kT2))
         * sin(ictx->kT * (ictx->xx - ictx->yx))
        );
    }
};

void I1(double x, double* real, double* imag) {
    if (x < 0) {
        // if argument is negative, switch it to positive and flip the sign of the imaginary part of the result
        x = -x;
        *imag = -(2 - cos(x) - gsl_sf_sinc(x * M_1_PI)) / x + 2 * gsl_sf_Si(x);
    }
    else {
        *imag = (2 - cos(x) - gsl_sf_sinc(x * M_1_PI)) / x - 2 * gsl_sf_Si(x);
    }
    *real = -2 * EULER_GAMMA - gsl_sf_sinc(x * M_1_PI) + (cos(x) - 1) / (x * x) + 2 * gsl_sf_Ci(x) - 2 * log(x);
}

class H14qq : public HardFactor {
public:
    term_type get_type() {
        return quadrupole;
    }
    double real_singular_contribution(IntegrationContext* ictx) {
        return -1/(4*M_PI*M_PI*4*M_PI*M_PI) * ictx->alphasbar * ictx->qqfactor / ictx->z2 * ictx->S4rst *
         // real part of the singular contribution
         4*M_PI * ictx->ctx->Nc * (1 + ictx->xi2) / ictx->xi
         * ((ictx->xx - ictx->bx)*(ictx->yx - ictx->bx) + (ictx->xy - ictx->by)*(ictx->yy - ictx->by))
            / ( ictx->s2 * ictx->t2 )
         * cos(ictx->kT * (ictx->xx / ictx->xi - ictx->yx - (1.0/ictx->xi - 1.0) * ictx->bx));
    }
    double imag_singular_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI*4*M_PI*M_PI) * ictx->alphasbar * ictx->qqfactor / ictx->z2 * ictx->S4rst *
         // imaginary part of the singular contribution
         4*M_PI * ictx->ctx->Nc * (1 + ictx->xi2) / ictx->xi
         * ((ictx->xx - ictx->bx)*(ictx->yx - ictx->bx) + (ictx->xy - ictx->by)*(ictx->yy - ictx->by))
            / ( ictx->s2 * ictx->t2 )
         * sin(ictx->kT * (ictx->xx / ictx->xi - ictx->yx - (1.0/ictx->xi - 1.0) * ictx->bx));
    }
    double real_delta_contribution(IntegrationContext* ictx) {
        double real, imag;
        I1(ictx->kT * (ictx->yx - ictx->bx), &real, &imag);
        return 1/(4*M_PI*M_PI*4*M_PI*M_PI) * ictx->alphasbar * ictx->qqfactor / ictx->z2 * (ictx->S4rst - ictx->ctx->gdist->S4(ictx->s2, ictx->s2, 0, ictx)) *
         // real part of the delta contribution
         4*M_PI * ictx->ctx->Nc * ictx->ctx->Sperp / ictx->t2
         * (cos(ictx->kT * (ictx->xx - ictx->yx)) * real + sin(ictx->kT * (ictx->xx - ictx->yx)) * imag);
    }
    double imag_delta_contribution(IntegrationContext* ictx) {
        double real, imag;
        I1(ictx->kT * (ictx->yx - ictx->bx), &real, &imag);
        return 1/(4*M_PI*M_PI*4*M_PI*M_PI) * ictx->alphasbar * ictx->qqfactor / ictx->z2 * (ictx->S4rst - ictx->ctx->gdist->S4(ictx->s2, ictx->s2, 0, ictx)) *
         // imaginary part of the delta contribution
         4*M_PI * ictx->ctx->Nc * ictx->ctx->Sperp / ictx->t2
         * (cos(ictx->kT * (ictx->xx - ictx->yx)) * imag - sin(ictx->kT * (ictx->xx - ictx->yx)) * real);
    }
};

class H14qqResidual : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    double real_delta_contribution(IntegrationContext* ictx) {
        return -1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->qqfactor / ictx->z2 * ictx->ctx->gdist->S4(ictx->r2, ictx->r2, 0, ictx) *
         // real part of the delta contribution
         ictx->ctx->Nc * ictx->ctx->Sperp * (2.5 - 2.0*M_PI*M_PI/3.0)
         * cos(ictx->kT * (ictx->xx - ictx->bx));
    }
    double imag_delta_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->qqfactor / ictx->z2 * ictx->ctx->gdist->S4(ictx->r2, ictx->r2, 0, ictx) *
         // imaginary part of the delta contribution
         ictx->ctx->Nc * ictx->ctx->Sperp * (2.5 - 2.0*M_PI*M_PI/3.0)
         * sin(ictx->kT * (ictx->xx - ictx->bx));
    }
};

class H02gg : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    double real_delta_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI) * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r *
         // real part of the hard factor
         cos(ictx->kT * (ictx->xx - ictx->yx)); // take angle of k_perp to be 0
    }
    double imag_delta_contribution(IntegrationContext* ictx) {
        return -1/(4*M_PI*M_PI) * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r *
         // imaginary part of the hard factor
         sin(ictx->kT * (ictx->xx - ictx->yx)); 
    }
};

class H12gg : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    double real_singular_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r *
         ictx->ctx->Nc * 2 * ictx->xi
           * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
           * (cos(ictx->kT * (ictx->xx - ictx->yx)) + cos(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi) / ictx->xi2);
    }
    double imag_singular_contribution(IntegrationContext* ictx) {
        return -1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r *
         ictx->ctx->Nc * 2 * ictx->xi
           * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
           * (sin(ictx->kT * (ictx->xx - ictx->yx)) + sin(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi) / ictx->xi2);
    }
    double real_normal_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r *
         ictx->ctx->Nc * 2 * (1.0 / ictx->xi - 1.0 + ictx->xi - ictx->xi2)
           * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
           * (cos(ictx->kT * (ictx->xx - ictx->yx)) + cos(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi) / ictx->xi2);
    }
    double imag_normal_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r *
         ictx->ctx->Nc * 2 * (1.0 / ictx->xi - 1.0 + ictx->xi - ictx->xi2)
           * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
           * (sin(ictx->kT * (ictx->xx - ictx->yx)) + sin(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi) / ictx->xi2);
    }
    double real_delta_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r * (
         ictx->ctx->Nc * (11.0/6.0 - 2 * ictx->ctx->Nf * ictx->ctx->TR / (3 * ictx->ctx->Nc))       // P_gg term
           * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
           * (cos(ictx->kT * (ictx->xx - ictx->yx)) + cos(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi) / ictx->xi2)
         - ictx->ctx->Nc * 2 * (11.0/6.0 - 2 * ictx->ctx->Nf * ictx->ctx->TR / (3 * ictx->ctx->Nc)) // independent term
           * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->kT2))
           * cos(ictx->kT * (ictx->xx - ictx->yx))
        );
    }
    double imag_delta_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r * (
         ictx->ctx->Nc * (11.0/6.0 - 2 * ictx->ctx->Nf * ictx->ctx->TR / (3 * ictx->ctx->Nc))       // P_gg term
           * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
           * (sin(ictx->kT * (ictx->xx - ictx->yx)) + sin(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi) / ictx->xi2)
         - ictx->ctx->Nc * 2 * (11.0/6.0 - 2 * ictx->ctx->Nf * ictx->ctx->TR / (3 * ictx->ctx->Nc)) // independent term
           * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->kT2))
           * sin(ictx->kT * (ictx->xx - ictx->yx))
        );
    }
};

void h12qqbar_internal_integrand(unsigned int ncoords, const double* coordinates, void* closure, unsigned int nvalues, double* values) {
    assert(ncoords == 1);
    assert(nvalues == 2);
    double xip = coordinates[0];
    double xip2 = xip * xip;
    double kr = *((double*)closure);
    values[0] = (2 * xip2 - 2 * xip + 1) * cos(kr);
    values[1] = -(2 * xip2 - 2 * xip + 1) * sin(kr);
}

void h12qqbar_internal_integral(double kr, double* real, double* imag) {
    double result[2];
    double error[2];
    double xmin1D[] = {0.0d};
    double xmax1D[] = {1.0d};
    int status = adapt_integrate(2, h12qqbar_internal_integrand, &kr, 1, xmin1D, xmax1D, 200000, 0, 1e-4, result, error);
    if (status != SUCCESS) {
        cerr << "Error in 1D integration (probably memory)" << endl;
        exit(1);
    }
    if (real) {
        *real = result[0];
    }
    if (imag) {
        *imag = result[1];
    }
}


// class H12qqbar : public HardFactor {
// public:
//     term_type get_type() {
//         return quadrupole;
//     }
//     double real_delta_contribution(IntegrationContext* ictx) {
//         return 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->ggfactor / ictx->z2 *
//          8 * M_PI * ictx->ctx->Nf * ictx->ctx->TR
//            * (ictx->ctx->gdist->S2(ictx->s2, ictx) * h12qqbar_internal_integral(ictx->kT * sqrt(ictx->r2)) - ictx->ctx->gdist->S2(ictx->t2, ictx) * h12qqbar_internal_integral(0)) * ictx->ctx->gdist->S2(ictx->t2, ictx)
//            / ictx->r2
//            * cos(ictx->kT * (ictx->yx - ictx->bx));
//     }
//     double imag_delta_contribution(IntegrationContext* ictx) {
//     }
// };

// class H12qqbarResidual : public HardFactor {
// public:
//     term_type get_type() {
//         return dipole;
//     }
//     double real_delta_contribution(IntegrationContext* ictx) {
//         return 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r *
//          8 * M_PI * ictx->ctx->Nf * ictx->ctx->TR
//            *// stuff...
//            * cos(ictx->kT * (ictx->yx - ictx->bx));
//     }
//     double imag_delta_contribution(IntegrationContext* ictx) {
//     }
// };

class H112gq : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    double real_normal_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->gqfactor / ictx->z2 * ictx->S2r *
         // real part of the normal contribution
         0.5 * ictx->ctx->Nc * (1 + (1 - ictx->xi)*(1 - ictx->xi)) / (ictx->xi * ictx->xi2) * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
         * cos(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi);
    }
    double imag_normal_contribution(IntegrationContext* ictx) {
        return -1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->gqfactor / ictx->z2 * ictx->S2r *
         // imaginary part of the normal contribution
         0.5 * ictx->ctx->Nc * (1 + (1 - ictx->xi)*(1 - ictx->xi)) / (ictx->xi * ictx->xi2) * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
         * sin(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi);
    }
};

class H122gq : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    double real_normal_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->gqfactor / ictx->z2 * ictx->S2r * ictx->S2r *
         // real part of the normal contribution
         0.5 * ictx->ctx->Nc * (1 + (1 - ictx->xi)*(1 - ictx->xi)) / ictx->xi * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
         * cos(ictx->kT * (ictx->xx - ictx->yx));
    }
    double imag_normal_contribution(IntegrationContext* ictx) {
        return -1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->gqfactor / ictx->z2 * ictx->S2r * ictx->S2r * // assumes S2 depends only on the magnitude of x - y
         // imaginary part of the normal contribution
         0.5 * ictx->ctx->Nc * (1 + (1 - ictx->xi)*(1 - ictx->xi)) / ictx->xi * (-2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
         * sin(ictx->kT * (ictx->xx - ictx->yx));
    }
};

class H14gq : public HardFactor {
public:
    term_type get_type() {
        return quadrupole;
    }
    double real_normal_contribution(IntegrationContext* ictx) {
        return -1/(4*M_PI*M_PI*4*M_PI*M_PI) * ictx->alphasbar * ictx->gqfactor / ictx->z2 * ictx->S4rst *
         // real part of the normal contribution
         4 * M_PI * ictx->ctx->Nc * (1 + (1 - ictx->xi)*(1 - ictx->xi)) / (ictx->xi2)
         * ((ictx->xx - ictx->yx)*(ictx->bx - ictx->yx) + (ictx->xy - ictx->yy)*(ictx->by - ictx->yy))
            / ( ictx->r2 * ictx->t2 )
         * cos(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi + ictx->kT * (ictx->yx - ictx->bx));
    }
    double imag_normal_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI*4*M_PI*M_PI) * ictx->alphasbar * ictx->gqfactor / ictx->z2 * ictx->S4rst *
         // imaginary part of the normal contribution
         4 * M_PI * ictx->ctx->Nc * (1 + (1 - ictx->xi)*(1 - ictx->xi)) / (ictx->xi2)
         * ((ictx->xx - ictx->yx)*(ictx->bx - ictx->yx) + (ictx->xy - ictx->yy)*(ictx->by - ictx->yy))
            / ( ictx->r2 * ictx->t2 )
         * sin(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi + ictx->kT * (ictx->yx - ictx->bx));
    }
};

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
    Integrator(Context* ctx, size_t hflen, HardFactor** hflist) : n_dipole_terms(0), n_quadrupole_terms(0), current_term_type(NONE), callback(NULL) {
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
    ~Integrator() {
        free(dipole_terms);
        dipole_terms = NULL;
        free(quadrupole_terms);
        quadrupole_terms = NULL;
        delete ictx;
    }
    void update(double z, double y, double rx, double ry) {
        // I need (z = z, y = y, xx = sx + bx, xy = sy + by, yx = tx + bx, yy = ty + by
        // we assume no impact parameter dependence, so I can arbitrarily set bx, by, tx, ty to 0
        assert(current_term_type == dipole);
        ictx->update(z, y, rx, ry, 0, 0, 0, 0);
    }
    void update(double z, double y, double sx, double sy, double tx, double ty) {
        // I need (z = z, y = y, xx = sx + bx, xy = sy + by, yx = tx + bx, yy = ty + by
        // we assume no impact parameter dependence, so I can arbitrarily set bx, by to 0
        assert(current_term_type == quadrupole);
        ictx->update(z, y, sx, sy, tx, ty, 0, 0);
    }
    void evaluate_1D_integrand(double* real, double* imag);
    void evaluate_2D_integrand(double* real, double* imag);
    void integrate(double* real, double* imag);
    void set_current_term_type(term_type new_term_type) {
        current_term_type = new_term_type;
    }
    void set_callback(void (*callback)(IntegrationContext*, double, double)) {
        this->callback = callback;
    }
};

void Integrator::evaluate_1D_integrand(double* real, double* imag) {
    double l_real = 0.0, l_imag = 0.0; // l for "local"
    double log_factor = log(1 - ictx->ctx->tau / ictx->z);
    assert(log_factor == log_factor);
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
        l_real += terms[i]->real_singular_contribution(ictx) * log_factor;
        assert(l_real == l_real);
        l_imag += terms[i]->imag_singular_contribution(ictx) * log_factor;
        assert(l_imag == l_imag);
        l_real += terms[i]->real_delta_contribution(ictx);
        assert(l_real == l_real);
        l_imag += terms[i]->imag_delta_contribution(ictx);
        assert(l_imag == l_imag);
    }
    if (callback) {
        callback(ictx, l_real, l_imag);
    }
    *real = l_real;
    *imag = l_imag;
}

void Integrator::evaluate_2D_integrand(double* real, double* imag) {
    double l_real = 0.0, l_imag = 0.0; // l for "local"
    double jacobian =  (1 - ictx->ctx->tau / ictx->z) / (1 - ictx->ctx->tau); // Jacobian from y to xi
    assert(jacobian == jacobian);
    double xi_factor = 1.0 / (1 - ictx->xi);
    assert(xi_factor == xi_factor);
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
        l_real += terms[i]->real_singular_contribution(ictx) * xi_factor;
        assert(l_real == l_real);
        l_imag += terms[i]->imag_singular_contribution(ictx) * xi_factor;
        assert(l_imag == l_imag);
        l_real += terms[i]->real_normal_contribution(ictx);
        assert(l_real == l_real);
        l_imag += terms[i]->imag_normal_contribution(ictx);
        assert(l_imag == l_imag);
    }
    if (callback) {
        callback(ictx, jacobian * l_real, jacobian * l_imag);
    }
    ictx->update(ictx->z, 1, ictx->xx, ictx->xy, ictx->yx, ictx->yy, ictx->bx, ictx->by);
    for (size_t i = 0; i < n_terms; i++) {
        l_real -= terms[i]->real_singular_contribution(ictx) * xi_factor;
        assert(l_real == l_real);
        l_imag -= terms[i]->imag_singular_contribution(ictx) * xi_factor;
        assert(l_imag == l_imag);
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
    gsl_monte_miser_integrate(&f, min, max, dim, 1000000, rng, s, p_result, p_abserr);
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
    if (callback) {
        (*callback)(p_result, p_abserr, s);
    }
    do {
        gsl_monte_vegas_integrate(&f, min, max, dim, 100000, rng, s, p_result, p_abserr);
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

void Integrator::integrate(double* real, double* imag) {
    double tmp_result;
    double tmp_error;
    double min1D[] = {ictx->ctx->tau, -inf, -inf, -inf, -inf};
    double max1D[] = {1.0, inf, inf, inf, inf};
    double min2D[] = {ictx->ctx->tau, ictx->ctx->tau, -inf, -inf, -inf, -inf};
    double max2D[] = {1.0, 1.0, inf, inf, inf, inf};
    double result = 0.0;
    double abserr = 0.0;
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
        abserr += tmp_error;
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
        abserr += tmp_error;
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
        abserr += tmp_error;
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
        abserr += tmp_error;
    }

    *real = result;
    *imag = 0;
}

#ifdef DATAGRID
void write_hf_values(Context* ctx) {
    double real, imag;
    GBWGluonDistribution* gdist = new GBWGluonDistribution();
    IntegrationContext* ictx = new IntegrationContext(ctx);
    H02qq* hf02qq = new H02qq();
    H12qq* hf12qq = new H12qq();
    for (double z = ctx->tau; z <= 1; z += (1 - ctx->tau) / 20) {
        for (double rx = -5; rx <= 5; rx += 20./270.) {
            for (double ry = -5; ry <= 5; ry += 20./120.) {
                ictx->update(z, 1, rx, ry, 0, 0, 0, 0);
                cout << z << "\t" << rx << "\t" << ry << "\t"
                     << hf02qq->real_delta_contribution(ictx) << "\t" << hf02qq->imag_delta_contribution(ictx) << "\t"
                     << hf12qq->real_delta_contribution(ictx) << "\t" << hf12qq->imag_delta_contribution(ictx) << endl;
            }
        }
    }
}

int main(int argc, char** argv) {
    Context gctx(
      0.000304, // x0
      197,      // A
      0.56,     // c
      0.288,    // lambda
      10,       // mu2
      3,        // Nc
      3,        // Nf
      1.5,      // CF
      1.0,      // Sperp
      0.4,      // pT2
      200,      // sqs
      3.2,      // Y
      0.2 / (2*M_PI), // alphasbar
      "mstw2008nlo.00.dat", "PINLO.DAT");
    write_hf_values(&gctx);
    return 0;
}
#else
void write_data_point(IntegrationContext* ictx, double real, double imag) {
    if (ictx) {
//         if ((++count) % 500 == 0) {
            cout
            << ictx->z << "\t"
            << ictx->xi << "\t"
            << ictx->xx << "\t"
            << ictx->xy << "\t"
            << ictx->yx << "\t"
            << ictx->yy << "\t"
            << ictx->bx << "\t"
            << ictx->by << "\t"
            << ictx->kT2 << "\t"
            << ictx->Qs2 << "\t"
            << ictx->xp << "\t"
            << ictx->xg << "\t"
            << ictx->qqfactor << "\t"
            << ictx->gqfactor << "\t"
            << real << "\t"
            << imag << endl;
//         }
    }
    else {
        cout << endl;
    }
}

double calculateHardFactor(Context* ctx, size_t hflen, HardFactor** hflist) {
    double real, imag;
    Integrator* integrator = new Integrator(ctx, hflen, hflist);
    if (trace) {
        integrator->set_callback(write_data_point);
    }
    integrator->integrate(&real, &imag);
    delete integrator;
    return real;
}

int main(int argc, char** argv) {
    GluonDistribution* gdist = new GBWGluonDistribution();
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--miser")==0) {
            integration_strategy = MC_MISER;
        }
        else if (strcmp(argv[i], "--vegas")==0) {
            integration_strategy = MC_VEGAS;
        }
        else if (strcmp(argv[i], "--trace")==0) {
            trace = true;
        }
    }
    gsl_rng_env_setup();
    Context gctx(
      0.000304, // x0
      197,      // A
      0.56,     // c
      0.288,    // lambda
      10,       // mu2
      3,        // Nc
      3,        // Nf
      1.5,      // CF
      0.5,      // TR
      1.0,      // Sperp
      1.0,      // pT2 (dummy value)
      200,      // sqs
      3.2,      // Y
      0.2 / (2*M_PI), // alphasbar
      gdist,
      "mstw2008nlo.00.dat", "PINLO.DAT");

//     double pT[] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4};
    double pT[] = {1.01};
    size_t pTlen = sizeof(pT)/sizeof(double);
    HardFactor* h02qq[] = {new H02qq()};
    HardFactor* h12qq[] = {new H12qq()};
    HardFactor* h14qq[] = {new H14qq(), new H14qqResidual()};
    HardFactor* h112gq[] = {new H112gq()};
    HardFactor* h122gq[] = {new H122gq()};
    HardFactor* h14gq[] = {new H14gq()};
    size_t hflen = 6;
    double yield[hflen*pTlen];
    
    for (size_t i = 0; i < pTlen; i++) {
        gctx.pT2 = pT[i]*pT[i];
        gctx.recalculate();
        cerr << "Beginning calculation at pT = " << pT[i] << endl;
        yield[hflen*i + 0] = calculateHardFactor(&gctx, 1, h02qq);
        yield[hflen*i + 1] = calculateHardFactor(&gctx, 1, h12qq);
        yield[hflen*i + 2] = calculateHardFactor(&gctx, 2, h14qq);
        yield[hflen*i + 3] = calculateHardFactor(&gctx, 1, h112gq);
        yield[hflen*i + 4] = calculateHardFactor(&gctx, 1, h122gq);
        yield[hflen*i + 5] = calculateHardFactor(&gctx, 1, h14gq);
        cerr << "...done" << endl;
    }
    cout << "pT\th02qq\th12qq\th14qq\th112gq\th122gq\th14gq" << endl;
    for (size_t i = 0; i < pTlen; i++) {
        cout << pT[i] << "\t";
        for (size_t j = 0; j < hflen; j++) {
            cout << yield[hflen*i + j] << "\t";
        }
        cout << endl;
    }
    delete gdist;
    return 0;
}
#endif
