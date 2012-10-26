#include <cassert>
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

class IntegrationContext;
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
    double Sperp;
    double pT2;
    double sqs;
    double Y;
    double Q02x0lambda;
    double tau;
    double alphasbar_fixed;

    // NOTE: these contain state that should be associated with IntegrationContext.
    // So don't use one Context with more than one IntegrationContext at once.
    c_mstwpdf* pdf_object;
    DSSpiNLO* ff_object;

    Context(double x0, double A, double c, double lambda, double mu2, double Nc, double Nf, double CF, double Sperp, double pT2, double sqs, double Y, double alphasbar_fixed, const char* pdf_filename, const char* ff_filename) :
     x0(x0), A(A), c(c), lambda(lambda), mu2(mu2), Nc(Nc), Nf(Nf), CF(CF), Sperp(Sperp), pT2(pT2), sqs(sqs), Y(Y), alphasbar_fixed(alphasbar_fixed) {
         Q02x0lambda = c * pow(A, 1.0d/3.0d) * pow(x0, lambda);
         tau = sqrt(pT2)/sqs*exp(Y);
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

    virtual double alphasbar(double kT2) {
        return alphasbar_fixed;
    }
};

class GluonDistribution {
public:
    virtual double S2(IntegrationContext* ictx) = 0;
    virtual double S4(IntegrationContext* ictx) = 0;
};

class GBWGluonDistribution: public GluonDistribution {
public:
    double S2(IntegrationContext* ictx);
    double S4(IntegrationContext* ictx);
};

// Not even anywhere close to thread-safe!
class IntegrationContext {
public:
    Context* ctx;
    GluonDistribution* gdist;
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
    double r2;
    double Qs2;
    double alphasbar;
    double quarkfactor;
    double S2, S4;
    
    IntegrationContext(Context* ctx, GluonDistribution* gdist) :
      ctx(ctx), gdist(gdist),
      z(0), xi(0),
      z2(0), xi2(0),
      xx(0), xy(0),
      yx(0), yy(0),
      bx(0), by(0),
      kT2(0), kT(0),
      xp(0), xg(0),
      Qs2(0), r2(0),
      alphasbar(0),
      quarkfactor(0),
      S2(0), S4(0) {
    };
    void update(double z, double y, double xx, double xy, double yx, double yy, double bx, double by);
};

void IntegrationContext::update(double z, double y, double xx, double xy, double yx, double yy, double bx, double by) {
    double real_singular = 0.0, imag_singular = 0.0,
           real_normal   = 0.0, imag_normal   = 0.0,
           real_delta    = 0.0, imag_delta    = 0.0;
    double quarkfactor = 0.0d;
    c_mstwpdf* pdf_object = ctx->pdf_object;
    DSSpiNLO* ff_object = ctx->ff_object;
    this->z = z;
    this->z2 = z*z;
    this->xi = (y * (z - ctx->tau) - ctx->tau * (z - 1)) / (z * (1 - ctx->tau));
    this->xi2 = xi*xi;
    this->xx = xx;
    this->xy = xy;
    this->yx = yx;
    this->yy = yy;
    this->bx = bx;
    this->by = by;
    this->r2 = (xx - yx) * (xx- yx) - (xy - yy) * (xy - yy);
    this->xp = ctx->tau / (z * xi);
    this->kT2 = ctx->pT2 / this->z2;
    this->kT = sqrt(this->kT2);
    this->xg = kT / ctx->sqs * exp(-ctx->Y);
    this->Qs2 = ctx->Q02x0lambda / pow(this->xg, ctx->lambda); // Q_0^2 (x_0 / x)^lambda
    this->alphasbar = ctx->alphasbar(this->kT2);

    // Calculate the new gluon distribution values
    // this has to be done after kinematics are updated
    this->S2 = gdist->S2(this);
    this->S4 = gdist->S4(this);

    // Calculate the new quark factor
    pdf_object->update(xp, sqrt(ctx->mu2));
    ff_object->update(z, ctx->mu2);
    
    // Proton contributions:
    // up quark
    quarkfactor += (pdf_object->cont.upv + pdf_object->cont.usea) * ff_object->fragmentation(DSSpiNLO::up, DSSpiNLO::pi_minus);
    // antiup quark
    quarkfactor += pdf_object->cont.usea * ff_object->fragmentation(DSSpiNLO::up_bar, DSSpiNLO::pi_minus);
    // down quark
    quarkfactor += (pdf_object->cont.dnv + pdf_object->cont.dsea) * ff_object->fragmentation(DSSpiNLO::down, DSSpiNLO::pi_minus);
    // antidown quark
    quarkfactor += pdf_object->cont.dsea * ff_object->fragmentation(DSSpiNLO::down_bar, DSSpiNLO::pi_minus);
    // strange quark
    quarkfactor += pdf_object->cont.str * ff_object->fragmentation(DSSpiNLO::strange, DSSpiNLO::pi_minus);
    // antistrange quark
    quarkfactor += pdf_object->cont.sbar * ff_object->fragmentation(DSSpiNLO::strange_bar, DSSpiNLO::pi_minus);
    
    // Neutron contributions (for deuteron collisions), assuming isospin symmetry:
    // up quark
    quarkfactor += (pdf_object->cont.dnv + pdf_object->cont.dsea) * ff_object->fragmentation(DSSpiNLO::up, DSSpiNLO::pi_minus);
    // antiup quark
    quarkfactor += pdf_object->cont.dsea * ff_object->fragmentation(DSSpiNLO::up_bar, DSSpiNLO::pi_minus);
    // down quark
    quarkfactor += (pdf_object->cont.upv + pdf_object->cont.usea) * ff_object->fragmentation(DSSpiNLO::down, DSSpiNLO::pi_minus);
    // antidown quark
    quarkfactor += pdf_object->cont.usea * ff_object->fragmentation(DSSpiNLO::down_bar, DSSpiNLO::pi_minus);
    // strange quark
    quarkfactor += pdf_object->cont.str * ff_object->fragmentation(DSSpiNLO::strange, DSSpiNLO::pi_minus);
    // antistrange quark
    quarkfactor += pdf_object->cont.sbar * ff_object->fragmentation(DSSpiNLO::strange_bar, DSSpiNLO::pi_minus);
    
    this->quarkfactor = quarkfactor;
}

double GBWGluonDistribution::S2(IntegrationContext* ictx) {
    return 1 / (M_PI * ictx->Qs2) * exp(-0.25 * ictx->r2 * ictx->Qs2);
}
double GBWGluonDistribution::S4(IntegrationContext* ictx) {
    return 1 / (M_PI * ictx->Qs2) * exp(-0.25 * (
        (ictx->xx - ictx->bx)*(ictx->xx - ictx->bx) + (ictx->xy - ictx->by)*(ictx->xy - ictx->by)
        + (ictx->yx - ictx->bx)*(ictx->yx - ictx->bx) + (ictx->yy - ictx->by)*(ictx->yy - ictx->by)
    ) * ictx->Qs2);
}

class HardFactor {
public:
    virtual double real_singular_contribution(IntegrationContext* ictx) { return 0; }
    virtual double imag_singular_contribution(IntegrationContext* ictx) { return 0; }
    virtual double real_normal_contribution(IntegrationContext* ictx) { return 0; }
    virtual double imag_normal_contribution(IntegrationContext* ictx) { return 0; }
    virtual double real_delta_contribution(IntegrationContext* ictx) { return 0; }
    virtual double imag_delta_contribution(IntegrationContext* ictx) { return 0; }
};

class H02qq : public HardFactor {
public:
    double real_delta_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI) * ictx->quarkfactor / ictx->z2 * ictx->S2 *
         // real part of the hard factor
         cos(ictx->kT * (ictx->xx - ictx->yx)); // take angle of k_perp to be 0
    }
    double imag_delta_contribution(IntegrationContext* ictx) {
        return -1/(4*M_PI*M_PI) * ictx->quarkfactor / ictx->z2 * ictx->S2 *
         // imaginary part of the hard factor
         sin(ictx->kT * (ictx->xx - ictx->yx)); 
    }
};

class H12qq : public HardFactor {
public:
    double real_singular_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->quarkfactor / ictx->z2 * ictx->S2 *
         // real part of the singular contribution
         ictx->ctx->CF * (1 + ictx->xi2) * (2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
         * (cos(ictx->kT * (ictx->xx - ictx->yx)) + cos(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi) / ictx->xi2);
    }
    double imag_singular_contribution(IntegrationContext* ictx) {
        return -1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->quarkfactor / ictx->z2 * ictx->S2 *
         // imaginary part of the singular contribution
         ictx->ctx->CF * (1 + ictx->xi2) * (2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
         * (sin(ictx->kT * (ictx->xx - ictx->yx)) + sin(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi) / ictx->xi2);
    }
    double real_delta_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->quarkfactor / ictx->z2 * ictx->S2 * (
         // real part of the delta contribution
         1.5 * ictx->ctx->CF * (2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
         * (cos(ictx->kT * (ictx->xx - ictx->yx)) + cos(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi) / ictx->xi2)
         - 3.0 * ictx->ctx->CF * (2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->kT2))
         * cos(ictx->kT * (ictx->xx - ictx->yx))
        );
    }
    double imag_delta_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI) * ictx->alphasbar * ictx->quarkfactor / ictx->z2 * ictx->S2 * (
         // imaginary part of the delta contribution
         1.5 * ictx->ctx->CF * (2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->ctx->mu2))
         * (sin(ictx->kT * (ictx->xx - ictx->yx)) + sin(ictx->kT * (ictx->xx - ictx->yx) / ictx->xi) / ictx->xi2)
         - 3.0 * ictx->ctx->CF * (2*EULER_GAMMA + log(4) - log(ictx->r2 * ictx->kT2))
         * sin(ictx->kT * (ictx->xx - ictx->yx))
        );
    }
};

class H14qq : public HardFactor {
public:
    double real_singular_contribution(IntegrationContext* ictx) {
        return -1/(4*M_PI*M_PI*4*M_PI*M_PI) * ictx->alphasbar * ictx->quarkfactor / ictx->z2 * ictx->S4 *
         // real part of the singular contribution
         4*M_PI * ictx->ctx->Nc * (1 + ictx->xi2) / ictx->xi
         * ((ictx->xx - ictx->bx)*(ictx->yx - ictx->bx) + (ictx->xy - ictx->by)*(ictx->yy - ictx->by))
            / ( ((ictx->xx - ictx->bx)*(ictx->xx - ictx->bx) + (ictx->xy - ictx->by)*(ictx->xy - ictx->by)) * ((ictx->yx - ictx->bx)*(ictx->yx - ictx->bx) + (ictx->yy - ictx->by)*(ictx->yy - ictx->by)) )
         * cos(ictx->kT * (ictx->xx / ictx->xi - ictx->yx - (1.0/ictx->xi - 1.0) * ictx->bx));
    }
    double imag_singular_contribution(IntegrationContext* ictx) {
        return 1/(4*M_PI*M_PI*4*M_PI*M_PI) * ictx->alphasbar * ictx->quarkfactor / ictx->z2 * ictx->S4 *
         // imaginary part of the singular contribution
         4*M_PI * ictx->ctx->Nc * (1 + ictx->xi2) / ictx->xi
         * ((ictx->xx - ictx->bx)*(ictx->yx - ictx->bx) + (ictx->xy - ictx->by)*(ictx->yy - ictx->by))
            / ( ((ictx->xx - ictx->bx)*(ictx->xx - ictx->bx) + (ictx->xy - ictx->by)*(ictx->xy - ictx->by)) * ((ictx->yx - ictx->bx)*(ictx->yx - ictx->bx) + (ictx->yy - ictx->by)*(ictx->yy - ictx->by)) )
         * sin(ictx->kT * (ictx->xx / ictx->xi - ictx->yx - (1.0/ictx->xi - 1.0) * ictx->bx));
    }
//     double real_delta_contribution(IntegrationContext* ictx) {
//         return -1/(4*M_PI*M_PI*4*M_PI*M_PI) * ictx->alphasbar * ictx->quarkfactor / ictx->z2 * ictx->S4 *
//          // real part of the delta contribution
//          1; // TODO: implement
//     }
};

class Integrator {
private:
    IntegrationContext* ictx;
    HardFactor** hflist;
    size_t hflen;
public:
    Integrator(Context* ctx, GluonDistribution* gdist, size_t hflen, HardFactor** hflist) : hflist(hflist), hflen(hflen) {
        ictx = new IntegrationContext(ctx, gdist);
    }
    ~Integrator() {
        delete ictx;
    }
    void update(double z, double y, double xx, double xy, double yx, double yy, double bx, double by) {
        ictx->update(z, y, xx, xy, yx, yy, bx, by);
    }
    void evaluate_1D_integrand(double* real, double* imag);
    void evaluate_2D_integrand(double* real, double* imag);
    void integrate(double* real, double* imag);
};

void Integrator::evaluate_1D_integrand(double* real, double* imag) {
    double l_real = 0.0, l_imag = 0.0; // l for "local"
    double log_factor = log(1 - ictx->ctx->tau / ictx->z);
    ictx->xi = 1;
    for (size_t i = 0; i < hflen; i++) {
        l_real += hflist[i]->real_singular_contribution(ictx) * log_factor;
        l_imag += hflist[i]->imag_singular_contribution(ictx) * log_factor;
        l_real += hflist[i]->real_delta_contribution(ictx);
        l_imag += hflist[i]->imag_delta_contribution(ictx);
    }
    *real = l_real;
    *imag = l_imag;
}

void Integrator::evaluate_2D_integrand(double* real, double* imag) {
    double l_real = 0.0, l_imag = 0.0; // l for "local"
    double jacobian =  (1 - ictx->ctx->tau / ictx->z) / (1 - ictx->ctx->tau); // Jacobian from y to xi
    double xi_factor = 1.0 / (1 - ictx->xi);
    for (size_t i = 0; i < hflen; i++) {
        l_real += hflist[i]->real_singular_contribution(ictx) * xi_factor;
        l_imag += hflist[i]->imag_singular_contribution(ictx) * xi_factor;
        l_real += hflist[i]->real_delta_contribution(ictx);
        l_imag += hflist[i]->imag_delta_contribution(ictx);
    }
    ictx->xi = 1;
    for (size_t i = 0; i < hflen; i++) {
        l_real -= hflist[i]->real_singular_contribution(ictx) * xi_factor;
        l_imag -= hflist[i]->imag_singular_contribution(ictx) * xi_factor;
    }
    *real = l_real;
    *imag = l_imag;
}

void cubature_wrapper_1D(unsigned int ncoords, const double* coordinates, void* closure, unsigned int nvalues, double* values) {
    Integrator* integrator = (Integrator*)closure;
    assert(ncoords == 7);
    assert(nvalues == 2);
    //                    z            y      xx              xy            yx               yy              bx              by
    integrator->update(coordinates[0], 1, coordinates[1], coordinates[2], coordinates[3], coordinates[4], coordinates[5], coordinates[6]);
    integrator->evaluate_2D_integrand(&(values[0]), &(values[1]));
}

void cubature_wrapper_2D(unsigned int ncoords, const double* coordinates, void* closure, unsigned int nvalues, double* values) {
    Integrator* integrator = (Integrator*)closure;
    assert(ncoords == 8);
    assert(nvalues == 2);
    //                    z               y               xx              xy              yx               yy              bx              by
    integrator->update(coordinates[0], coordinates[1], coordinates[2], coordinates[3], coordinates[4], coordinates[5], coordinates[6], coordinates[7]);
    integrator->evaluate_2D_integrand(&(values[0]), &(values[1]));
}


void Integrator::integrate(double* real, double* imag) {
    double l_real = 0.0, l_imag = 0.0;
    double tmp_result[] = {0.0, 0.0};
    double tmp_error[] = {0.0, 0.0};
    double min1D[] = {ictx->ctx->tau, -100, -100, -100, -100, -100, -100};
    double max1D[] = {1.0, 100, 100, 100, 100, 100, 100};
    double min2D[] = {ictx->ctx->tau, ictx->ctx->tau, -100, -100, -100, -100, -100, -100};
    double max2D[] = {1.0, 1.0, 100, 100, 100, 100, 100, 100};
    int status = 0;
    // the 2D integration
    status = adapt_integrate(2, cubature_wrapper_2D, this, 8, min2D, max2D, 100000, 0, 1e-2, tmp_result, tmp_error);
    if (status != SUCCESS) {
        cerr << "Error in 2D integration (probably memory)" << endl;
        exit(1);
    }
    l_real += tmp_result[0];
    l_imag += tmp_result[1];
    // the 1D integration
    status = adapt_integrate(2, cubature_wrapper_1D, this, 7, min1D, max1D, 100000, 0, 1e-2, tmp_result, tmp_error);
    if (status != SUCCESS) {
        cerr << "Error in 1D integration (probably memory)" << endl;
        exit(1);
    }
    l_real += tmp_result[0];
    l_imag += tmp_result[1];
    *real = l_real;
    *imag = l_imag;
}

double calculateLOterm(Context* ctx) {
    double real, imag;
    GBWGluonDistribution* gdist = new GBWGluonDistribution();
    HardFactor* hflist[1];
    size_t hflen = sizeof(hflist)/sizeof(hflist[0]);
    Integrator* integrator = new Integrator(ctx, gdist, hflen, hflist);
    hflist[0] = new H02qq();
    integrator->integrate(&real, &imag);
    delete integrator;
    for (size_t i = 0; i < hflen; i++) {
        delete hflist[i];
    }
    return real;
}

double calculateNLOterm(Context* ctx) {
    // this only calculates part of the NLO term
    double real, imag;
    GBWGluonDistribution* gdist = new GBWGluonDistribution();
    HardFactor* hflist[1];
    size_t hflen = sizeof(hflist)/sizeof(hflist[0]);
    Integrator* integrator = new Integrator(ctx, gdist, hflen, hflist);
    hflist[0] = new H12qq();
    integrator->integrate(&real, &imag);
    delete integrator;
    for (size_t i = 0; i < hflen; i++) {
        delete hflist[i];
    }
    return real;
}

void fillYieldArray(double sqs, double Y, int pTlen, double* pT, double* yield) {
    int i;
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
      1.0,      // pT2 (dummy value)
      sqs,
      Y,
      0.2 / (2*M_PI), // alphasbar
      "mstw2008nlo.00.dat", "PINLO.DAT");
    for (i = 0; i < pTlen; i++) {
        gctx.pT2 = pT[i]*pT[i];
        cerr << "Beginning calculation at pT = " << pT[i] << endl;
        yield[2*i] = calculateLOterm(&gctx);
        yield[2*i+1] = calculateNLOterm(&gctx);
        cerr << "...done" << endl;
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

