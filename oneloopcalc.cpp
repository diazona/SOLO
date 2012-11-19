/*
 * A calculation of the NLO cross section of pA->pion collisions
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
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_sys.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include "cubature.h"
#include "mstwpdf.h"
#include "dss_pinlo.h"

#define checkfinite(d) assert(gsl_finite(d))

using namespace std;

const int SUCCESS = 0;

static enum {MC_PLAIN, MC_MISER, MC_VEGAS} integration_strategy = MC_VEGAS;
typedef enum {NONE=0, dipole=2, quadrupole=4} term_type;
static bool trace = false;

void eprint(char* s) { cerr << s << endl; }

class IntegrationContext;
class GluonDistribution;
class Coupling;
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
    
    GluonDistribution* gdist;
    Coupling* cpl;

    // NOTE: these contain state that should be associated with IntegrationContext.
    // So don't use one Context with more than one IntegrationContext at once.
    c_mstwpdf* pdf_object;
    DSSpiNLO* ff_object;

    Context(double x0, double A, double c, double lambda, double mu2, double Nc, double Nf, double CF, double TR, double Sperp, double pT2, double sqs, double Y, GluonDistribution* gdist, Coupling* cpl, const char* pdf_filename, const char* ff_filename) :
     x0(x0), A(A), c(c), lambda(lambda), mu2(mu2), Nc(Nc), Nf(Nf), CF(CF), TR(TR), Sperp(Sperp), pT2(pT2), sqs(sqs), Y(Y), gdist(gdist), cpl(cpl) {
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
};

class Coupling {
public:
    virtual double alphasbar(double kT2) = 0;
};

class FixedCoupling : public Coupling {
private:
    double value;
public:
    FixedCoupling(double alphasbar) : value(alphasbar) {}
    double alphasbar(double kT2) {
        return value;
    }
};

class LORunningCoupling : public Coupling {
private:
    double log_LambdaQCD;
    double inverse_beta_2pi;
    double regulator; // position of the Landau pole
public:
    LORunningCoupling(double LambdaQCD, double beta, double regulator) : log_LambdaQCD(log(LambdaQCD)), inverse_beta_2pi(0.5 / (M_PI * beta)), regulator(regulator) {}
    double alphasbar(double kT2) {
        return inverse_beta_2pi / (log(kT2 + regulator) - log_LambdaQCD);
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
    double rx, ry, r2;
    double sx, sy, s2;
    double tx, ty, t2;
    double Qs2;
    double alphasbar;
    double qqfactor;
    double ggfactor;
    double gqfactor;
    double qgfactor;
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
      rx(0), ry(0), r2(0),
      sx(0), sy(0), s2(0),
      tx(0), ty(0), t2(0),
      Qs2(0),
      alphasbar(0),
      qqfactor(0),
      ggfactor(0),
      gqfactor(0),
      qgfactor(0),
      S2r(0), S4rst(0) {
    };
    void update(double z, double y, double xx, double xy, double yx, double yy, double bx, double by);
};

void IntegrationContext::update(double z, double y, double xx, double xy, double yx, double yy, double bx, double by) {
    double qqfactor = 0.0d, ggfactor = 0.0d, gqfactor = 0.0d, qgfactor = 0.0d;
    c_mstwpdf* pdf_object = ctx->pdf_object;
    DSSpiNLO* ff_object = ctx->ff_object;
    assert(z <= 1);
    assert(z >= ctx->tau);
    assert(y <= 1);
    assert(y >= ctx->tau);
    checkfinite(xx);
    checkfinite(xy);
    checkfinite(yx);
    checkfinite(yy);
    checkfinite(bx);
    checkfinite(by);
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
    this->rx = this->xx - this->yx;
    this->ry = this->xy - this->yy;
    this->r2 = rx*rx + ry*ry;
    this->sx = this->xx - this->bx;
    this->sy = this->xy - this->by;
    this->s2 = sx*sx + sy*sy;
    this->tx = this->yx - this->bx;
    this->ty = this->yy - this->by;
    this->t2 = tx*tx + ty*ty;
    this->xp = ctx->tau / (z * xi);
    this->kT2 = ctx->pT2 / this->z2;
    this->kT = sqrt(this->kT2);
    this->xg = kT / ctx->sqs * exp(-ctx->Y);
    this->Qs2 = ctx->Q02x0lambda / pow(this->xg, ctx->lambda); // Q_0^2 (x_0 / x)^lambda
    this->alphasbar = ctx->cpl->alphasbar(this->kT2);

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

    
    // Proton contributions:
    qgfactor += pdf_object->cont.glu * (  ff_object->fragmentation(DSSpiNLO::up, DSSpiNLO::pi_minus)
                                        + ff_object->fragmentation(DSSpiNLO::up_bar, DSSpiNLO::pi_minus)
                                        + ff_object->fragmentation(DSSpiNLO::down, DSSpiNLO::pi_minus)
                                        + ff_object->fragmentation(DSSpiNLO::down_bar, DSSpiNLO::pi_minus)
                                        + ff_object->fragmentation(DSSpiNLO::strange, DSSpiNLO::pi_minus)
                                        + ff_object->fragmentation(DSSpiNLO::strange_bar, DSSpiNLO::pi_minus));
    
    // Neutron contributions (for deuteron collisions), assuming isospin symmetry:
    qgfactor *= 2;
    
    this->qgfactor = qgfactor;

    
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
    virtual void Fs(IntegrationContext* ictx, double* real, double* imag) { *real = 0; *imag = 0; }
    virtual void Fn(IntegrationContext* ictx, double* real, double* imag) { *real = 0; *imag = 0; }
    virtual void Fd(IntegrationContext* ictx, double* real, double* imag) { *real = 0; *imag = 0; }
};

class H02qq : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = 1/(4*M_PI*M_PI) * ictx->qqfactor / ictx->z2 * ictx->S2r;
        double phase = -ictx->kT * (ictx->xx - ictx->yx); // take angle of k_perp to be 0
        *real = amplitude * cos(phase);
        *imag = amplitude * sin(phase);
    }
};

class H12qq : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fs(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->qqfactor / ictx->z2 * ictx->S2r *
          ictx->ctx->CF * (1 + ictx->xi2) * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->ctx->mu2));
        double phase1 = -ictx->kT * (ictx->xx - ictx->yx);
        double phase2 = -ictx->kT * (ictx->xx - ictx->yx) / ictx->xi;
        *real = amplitude * (cos(phase1) + cos(phase2) / ictx->xi2);
        *imag = amplitude * (sin(phase1) + sin(phase2) / ictx->xi2);
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->qqfactor / ictx->z2 * ictx->S2r;
        double term1 = 1.5 * ictx->ctx->CF * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->ctx->mu2));
        double term2 = -3.0 * ictx->ctx->CF * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->kT2));
        double phase1 = ictx->kT * (ictx->xx - ictx->yx);
        double phase2 = ictx->kT * (ictx->xx - ictx->yx) / ictx->xi;
        *real = amplitude * (term1 * (cos(phase1) + cos(phase2) / ictx->xi2) + term2 * cos(phase1));
        *imag = amplitude * (term1 * (sin(phase1) + sin(phase2) / ictx->xi2) + term2 * sin(phase1));
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
    *real = -2 * M_EULER - gsl_sf_sinc(x * M_1_PI) + (cos(x) - 1) / (x * x) + 2 * gsl_sf_Ci(x) - 2 * log(x);
}

class H14qq : public HardFactor {
public:
    term_type get_type() {
        return quadrupole;
    }
    void Fs(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = -1/(4*M_PI*M_PI*4*M_PI*M_PI) * ictx->alphasbar * ictx->qqfactor / ictx->z2 * ictx->S4rst *
         4*M_PI * ictx->ctx->Nc * (1 + ictx->xi2) / ictx->xi
         * ((ictx->xx - ictx->bx)*(ictx->yx - ictx->bx) + (ictx->xy - ictx->by)*(ictx->yy - ictx->by))
            / ( ictx->s2 * ictx->t2 );
        double phase = -ictx->kT * (ictx->xx / ictx->xi - ictx->yx - (1.0/ictx->xi - 1.0) * ictx->bx);
        *real = amplitude * cos(phase);
        *imag = amplitude * sin(phase);
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag) {
        double realI, imagI;
        I1(ictx->kT * (ictx->yx - ictx->bx), &realI, &imagI);
        double amplitude = 1/(4*M_PI*M_PI*4*M_PI*M_PI) * ictx->alphasbar * ictx->qqfactor / ictx->z2 * (ictx->S4rst - ictx->ctx->gdist->S4(ictx->s2, ictx->s2, 0, ictx)) *
         4*M_PI * ictx->ctx->Nc * ictx->ctx->Sperp / ictx->t2;
        double phase = ictx->kT * (ictx->xx - ictx->yx);
        *real = amplitude * (cos(phase) * realI + sin(phase) * imagI);
        *imag = amplitude * (cos(phase) * imagI - sin(phase) * realI);
    }
};

class H14qqResidual : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = -1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->qqfactor / ictx->z2 * ictx->ctx->gdist->S4(ictx->r2, ictx->r2, 0, ictx) *
         ictx->ctx->Nc * ictx->ctx->Sperp * (2.5 - 2.0*M_PI*M_PI/3.0);
        double phase = ictx->kT * (ictx->xx - ictx->bx);
        *real = amplitude * cos(phase);
        *imag = amplitude * sin(phase);
    }
};

class H02gg : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = 1/(4*M_PI*M_PI) * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r;
        double phase = -ictx->kT * (ictx->xx - ictx->yx);
        *real = amplitude * cos(phase);
        *imag = amplitude * sin(phase);
    }
};

class H12gg : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fs(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r *
         ictx->ctx->Nc * 2 * ictx->xi
           * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->ctx->mu2));
        double phase1 = -ictx->kT * (ictx->xx - ictx->yx);
        double phase2 = -ictx->kT * (ictx->xx - ictx->yx) / ictx->xi;
        *real = amplitude * (cos(phase1) + cos(phase2) / ictx->xi2);
        *imag = amplitude * (sin(phase1) + sin(phase2) / ictx->xi2);
    }
    void Fn(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r *
         ictx->ctx->Nc * 2 * (1.0 / ictx->xi - 1.0 + ictx->xi - ictx->xi2)
           * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->ctx->mu2));
        double phase1 = ictx->kT * (ictx->xx - ictx->yx);
        double phase2 = ictx->kT * (ictx->xx - ictx->yx) / ictx->xi;
        *real = amplitude * (cos(phase1) + cos(phase2) / ictx->xi2);
        *imag = amplitude * (sin(phase1) + sin(phase2) / ictx->xi2);
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r;
        double term1 = ictx->ctx->Nc * (11.0/6.0 - 2 * ictx->ctx->Nf * ictx->ctx->TR / (3 * ictx->ctx->Nc))       // P_gg term
           * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->ctx->mu2));
        double term2 = -ictx->ctx->Nc * 2 * (11.0/6.0 - 2 * ictx->ctx->Nf * ictx->ctx->TR / (3 * ictx->ctx->Nc)) // independent term
           * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->kT2));
        double phase1 = ictx->kT * (ictx->xx - ictx->yx);
        double phase2 = ictx->kT * (ictx->xx - ictx->yx) / ictx->xi;
        *real = amplitude * (term1 * (cos(phase1) + cos(phase2) / ictx->xi2) + term2 * cos(phase1));
        *imag = amplitude * (term1 * (sin(phase1) + sin(phase2) / ictx->xi2) + term2 * cos(phase1));
    }
};

void I2(double x, double* real, double* imag) {
    double sx = sin(x), cx = cos(x);
    double ix = 1.0 / x;
    double ix2 = ix * ix;
    double ix3 = ix2 * ix;
    *real = -4 * ix3 * sx + 2 * ix2 * (1 + cx) + ix * sx;
    *imag = 4 * ix3 * (1 - cx) - 2 * ix2 * sx - ix * (1 - cx);
}

class H12qqbar : public HardFactor {
public:
    term_type get_type() {
        return quadrupole;
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag) {
        double realI, imagI;
        I2(ictx->kT * sqrt(ictx->r2), &realI, &imagI);
        double amplitude = 1/(4*M_PI*M_PI*4*M_PI*M_PI) * ictx->alphasbar * ictx->ggfactor / ictx->z2 *
         8 * M_PI * ictx->ctx->Nf * ictx->ctx->TR * ictx->ctx->Sperp
           * (ictx->ctx->gdist->S2(ictx->s2, ictx) - ictx->ctx->gdist->S2(ictx->t2, ictx)) * ictx->ctx->gdist->S2(ictx->t2, ictx)
           / ictx->r2;
        double phase = ictx->kT * (ictx->yx - ictx->bx);
        *real = amplitude * (cos(phase) * realI + sin(phase) * imagI);
        *imag = amplitude * (cos(phase) * imagI - sin(phase) * realI);
    }
};

class H12qqbarResidual : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->ggfactor / ictx->z2 * ictx->ctx->Sperp * ictx->S2r * ictx->S2r *
         4 * M_PI * ictx->ctx->Nf * ictx->ctx->TR * (13.0 / 18.0);
        double phase = -ictx->kT * (ictx->yx - ictx->bx);
        *real = amplitude * cos(phase);
        *imag = amplitude * sin(phase);
    }
};

void I3(double x, double* real, double* imag) {
    bool negative = false;
    if (x == 0) {
        *real = -11.0 / 12.0;
        *imag = 0;
        return;
    }
    else if (x < 0) {
        x = -x;
        negative = true;
    }
    double sx = sin(x), cx = cos(x);
    double si = gsl_sf_Si(x), ci = gsl_sf_Ci(x);
    double ix = 1.0 / x;
    double ix2 = ix * ix;
    double ix3 = ix2 * ix;
    *real = -0.5 * (1 + cx) * ix2 + cx * (ci - log(x) - M_EULER) + sx * (si - ix + ix3);
    *imag = (negative ? -1 : 1) * ( (1 - cx) * (ix - ix3) + cx * si + sx * (log(x) + 0.5 * ix2 + M_EULER - ci) );
}

class H16gg : public HardFactor {
public:
    term_type get_type() {
        return quadrupole;
    }
    void Fs(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = M_1_PI*M_1_PI*M_1_PI * ictx->alphasbar * ictx->ctx->Nc * ictx->ctx->Sperp / ictx->z2
          * ictx->ggfactor * ictx->S4rst * (1 - ictx->xi + ictx->xi2)*(1 - ictx->xi + ictx->xi2) / ictx->xi2
          * ((ictx->xx - ictx->yx)*(ictx->yx - ictx->bx) + (ictx->xy - ictx->yy)*(ictx->yy - ictx->by)) / (ictx->r2 * ictx->t2);
        double phase = -ictx->kT * (ictx->xx - ictx->yx + (ictx->yx - ictx->bx) / ictx->xi);
        *real = amplitude * cos(phase);
        *imag = amplitude * sin(phase);
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag) {
        double realI, imagI;
        I3(ictx->kT * (ictx->yx - ictx->bx), &realI, &imagI);
        double amplitude = M_1_PI*M_1_PI*M_1_PI * ictx->alphasbar * ictx->ctx->Nc * ictx->ctx->Sperp / ictx->z2
          * ictx->ggfactor * ictx->ctx->gdist->S2(ictx->s2, ictx) * (ictx->ctx->gdist->S2(ictx->t2, ictx) * ictx->ctx->gdist->S2(ictx->r2, ictx) - ictx->ctx->gdist->S2(ictx->s2, ictx) * ictx->ctx->gdist->S2(0, ictx))
          / ictx->t2;
        double phase = -ictx->kT * (ictx->xx - ictx->yx);
        *real = amplitude * (cos(phase) * realI - sin(phase) * imagI);
        *imag = amplitude * (cos(phase) * imagI + sin(phase) * realI);
    }
};

class H16ggResidual : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = (-67.0 / 36.0 * M_1_PI*M_1_PI + 1.0 / 3.0) * ictx->alphasbar * ictx->ctx->Nc * ictx->ctx->Sperp / ictx->z2
          * ictx->ggfactor * ictx->S2r * ictx->S2r * ictx->ctx->gdist->S2(0, ictx);
        double phase = -ictx->kT * (ictx->xx - ictx->yx);
        *real = amplitude * cos(phase);
        *imag = amplitude * cos(phase);
    }
};

class H112gq : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fn(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->gqfactor / ictx->z2 * ictx->S2r *
         0.5 * ictx->ctx->Nc * (1 + (1 - ictx->xi)*(1 - ictx->xi)) / (ictx->xi * ictx->xi2) * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->ctx->mu2));
        double phase = -ictx->kT * (ictx->xx - ictx->yx) / ictx->xi;
        checkfinite(amplitude);
        checkfinite(phase);
        *real = amplitude * cos(phase);
        *imag = amplitude * sin(phase);
    }
};

class H122gq : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fn(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->gqfactor / ictx->z2 * ictx->S2r * ictx->S2r *
         0.5 * ictx->ctx->Nc * (1 + (1 - ictx->xi)*(1 - ictx->xi)) / ictx->xi * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->ctx->mu2));
        double phase = -ictx->kT * (ictx->xx - ictx->yx);
        checkfinite(amplitude);
        checkfinite(phase);
        *real = amplitude * cos(phase);
        *imag = amplitude * sin(phase);
    }
};

class H14gq : public HardFactor {
public:
    term_type get_type() {
        return quadrupole;
    }
    void Fn(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = -1/(4*M_PI*M_PI*4*M_PI*M_PI) * ictx->alphasbar * ictx->gqfactor / ictx->z2 * ictx->S4rst *
         4 * M_PI * ictx->ctx->Nc * (1 + (1 - ictx->xi)*(1 - ictx->xi)) / (ictx->xi2)
         * ((ictx->xx - ictx->yx)*(ictx->bx - ictx->yx) + (ictx->xy - ictx->yy)*(ictx->by - ictx->yy))
            / ( ictx->r2 * ictx->t2 );
        double phase = -(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi + ictx->kT * (ictx->yx - ictx->bx));
        checkfinite(amplitude);
        checkfinite(phase);
        *real = amplitude * cos(phase);
        *imag = amplitude * sin(phase);
    }
};

class H112qg : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fn(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = 0.125*M_1_PI*M_1_PI * ictx->alphasbar * ictx->ctx->Sperp / ictx->z2 * ictx->qgfactor * ictx->S2r *
         (1 - 2 * ictx->xi + 2 * ictx->xi2) * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->ctx->mu2));
        double phase = -ictx->kT * (ictx->xx - ictx->yx);
        checkfinite(amplitude);
        checkfinite(phase);
        *real = amplitude * cos(phase);
        *imag = amplitude * sin(phase);
    }
};

class H122qg : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fn(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = 0.125*M_1_PI*M_1_PI * ictx->alphasbar * ictx->ctx->Sperp / ictx->z2 * ictx->qgfactor * ictx->S2r *
         (1 / ictx->xi2 - 2 / ictx->xi + 2) * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->ctx->mu2));
        double phase = -ictx->kT * (ictx->xx - ictx->yx) / ictx->xi;
        checkfinite(amplitude);
        checkfinite(phase);
        *real = amplitude * cos(phase);
        *imag = amplitude * sin(phase);
    }
};

class H14qg : public HardFactor {
public:
    term_type get_type() {
        return quadrupole;
    }
    void Fn(IntegrationContext* ictx, double* real, double* imag) {
        double amplitude = 0.25*M_1_PI*M_1_PI*M_1_PI * ictx->alphasbar * ictx->ctx->Sperp / ictx->z2 * ictx->qgfactor * ictx->S4rst *
         (1 / ictx->xi - 2 + 2 * ictx->xi)
         * ((ictx->xx - ictx->yx)*(ictx->yx - ictx->bx) + (ictx->xy - ictx->yy)*(ictx->yy - ictx->by))
            / ( ictx->r2 * ictx->t2 );
        double phase = -(ictx->kT * (ictx->xx - ictx->yx) + ictx->kT * (ictx->yx - ictx->bx) / ictx->xi);
        checkfinite(amplitude);
        checkfinite(phase);
        *real = amplitude * cos(phase);
        *imag = amplitude * sin(phase);
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
    bool separate = false;
    GluonDistribution* gdist = new GBWGluonDistribution();
    Coupling* cpl = new FixedCoupling(0.2 / (2*M_PI));
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
        else if (strcmp(argv[i], "--separate")==0) {
            separate = true;
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
      gdist,
      cpl,
      "mstw2008nlo.00.dat", "PINLO.DAT");

    double pT[] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4};
    size_t pTlen = sizeof(pT)/sizeof(double);
    HardFactor* hflist[] = {
        new H02qq(), // 0
        new H02gg(), // 1
        new H12qq(), // 2
        new H14qq(), new H14qqResidual(), // 3,4
        new H12gg(), // 5
        new H12qqbar(), new H12qqbarResidual(), // 6,7
        new H16gg(), new H16ggResidual(), // 8,9
        new H112gq(), // 10
        new H122gq(), // 11
        new H14gq(), // 12
        new H112qg(), // 13
        new H122qg(), // 14
        new H14qg()}; // 15
    size_t hflen = 13;
    double yield[hflen*pTlen];
    
    for (size_t i = 0; i < pTlen; i++) {
        gctx.pT2 = pT[i]*pT[i];
        gctx.recalculate();
        cerr << "Beginning calculation at pT = " << pT[i] << endl;
        if (separate) {
            yield[hflen*i + 0] = calculateHardFactor(&gctx, 1, hflist + 0);
            yield[hflen*i + 1] = calculateHardFactor(&gctx, 1, hflist + 1);
            yield[hflen*i + 2] = calculateHardFactor(&gctx, 1, hflist + 2);
            yield[hflen*i + 3] = calculateHardFactor(&gctx, 2, hflist + 3);
            yield[hflen*i + 4] = calculateHardFactor(&gctx, 1, hflist + 5);
            yield[hflen*i + 5] = calculateHardFactor(&gctx, 2, hflist + 6);
            yield[hflen*i + 6] = calculateHardFactor(&gctx, 2, hflist + 8);
            yield[hflen*i + 7] = calculateHardFactor(&gctx, 1, hflist + 10);
            yield[hflen*i + 8] = calculateHardFactor(&gctx, 1, hflist + 11);
            yield[hflen*i + 9] = calculateHardFactor(&gctx, 1, hflist + 12);
            yield[hflen*i + 10] = calculateHardFactor(&gctx, 1, hflist + 13);
            yield[hflen*i + 11] = calculateHardFactor(&gctx, 1, hflist + 14);
            yield[hflen*i + 12] = calculateHardFactor(&gctx, 1, hflist + 15);
        }
        else {
            yield[2*i + 0] = calculateHardFactor(&gctx, 2, hflist); // leading order
            yield[2*i + 1] = calculateHardFactor(&gctx, 14, hflist+2); // next-to-leading order
        }
        cerr << "...done" << endl;
    }
    if (separate) {
        cout << "pT\th02qq\th02gg\th12qq\th14qq\th12gg\th12qqbar\th16gg\th112gq\th122gq\th14gq\th112qg\th122qg\th14qg\tlo\tnlo\tlo+nlo" << endl;
    }
    else {
        cout << "pT\tlo\tnlo\tlo+nlo" << endl;
    }
    for (size_t i = 0; i < pTlen; i++) {
        if (separate) {
            cout << pT[i] << "\t";

            double lo = 0, nlo = 0;
            size_t j;
            for (j = 0; j < 2; j++) {
                cout << yield[hflen*i + j] << "\t";
                lo += yield[hflen*i + j];
            }
            for (; j < hflen; j++) {
                cout << yield[hflen*i + j] << "\t";
                nlo += yield[hflen*i + j];
            }
            cout << lo << "\t" << nlo << "\t" << lo+nlo << endl;
        }
        else {
            cout << pT[i] << "\t" << yield[2*i+0] << "\t" << yield[2*i+1] << "\t" << yield[2*i+0]+yield[2*i+1] << endl;
        }
    }
    delete gdist;
    return 0;
}
#endif
