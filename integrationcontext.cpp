#include <cassert>
#include <gsl/gsl_sys.h>
#include "mstwpdf.h"
#include "dss_pinlo.h"
#include "integrationcontext.h"
#include "gluondist.h"

#define checkfinite(d) assert(gsl_finite(d))

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
    this->z2 = this->z*this->z;
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
    this->xp = ctx->tau / (this->z * xi);
    this->kT2 = ctx->pT2 / this->z2;
    this->kT = sqrt(this->kT2);
    this->xg = kT / ctx->sqs * exp(-ctx->Y);
    this->Qs2 = ctx->Q02x0lambda / pow(this->xg, ctx->lambda); // Q_0^2 (x_0 / x)^lambda
    this->alphasbar = ctx->cpl->alphasbar(this->kT2);

    // Calculate the new gluon distribution values
    // this has to be done after kinematics are updated
    this->S2r = ctx->gdist->S2(r2, this->Qs2);
    this->S4rst = ctx->gdist->S4(r2, s2, t2, this->Qs2);

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
