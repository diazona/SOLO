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
#include <gsl/gsl_sys.h>
#include "mstwpdf.h"
#include "dss_pinlo.h"
#include "integrationcontext.h"
#include "gluondist.h"

#define checkfinite(d) assert(gsl_finite(d))

void IntegrationContext::update_parton_factors(double z, double y) {
    assert(z <= 1);
    assert(z >= ctx->tau);
    assert(y <= 1);
    assert(y >= ctx->tau);
    this->z = z;
    this->z2 = this->z*this->z;
    if (y == 1.0d) {
        this->xi = 1.0d; // avoid floating-point roundoff error
    }
    else {
        this->xi = (y * (z - ctx->tau) - ctx->tau * (z - 1)) / (z * (1 - ctx->tau));
    }
    this->xi2 = xi*xi;
    this->xp = ctx->tau / (this->z * xi);
    this->kT2 = ctx->pT2 / this->z2;
    this->kT = sqrt(this->kT2);
    this->xg = kT / ctx->sqs * exp(-ctx->Y);
    const SaturationScale& satscale = ctx->gdist->satscale;
    this->Yg = satscale.Yx(this->xg);
    this->Qs2 = satscale.Qs2x(this->xg);
    this->Fk = ctx->gdist->F(this->kT2, this->Yg);
    this->alphas = ctx->cpl->alphas(this->kT2);
    this->alphas_2pi = this->alphas * 0.5 * M_1_PI;

    double qqfactor = 0.0d, ggfactor = 0.0d, gqfactor = 0.0d, qgfactor = 0.0d;
    c_mstwpdf* pdf_object = tlctx->pdf_object;
    DSSpiNLO* ff_object = tlctx->ff_object;
    
    DSSpiNLO::hadron hadron = ctx->hadron;
    
    // Calculate the new quark/gluon factors
    pdf_object->update(xp, sqrt(ctx->mu2));
    ff_object->update(z, ctx->mu2);
    
    // Proton contributions:
    qqfactor += (pdf_object->cont.upv + pdf_object->cont.usea) * ff_object->fragmentation(DSSpiNLO::up, hadron);
    qqfactor += pdf_object->cont.usea * ff_object->fragmentation(DSSpiNLO::up_bar, hadron);
    qqfactor += (pdf_object->cont.dnv + pdf_object->cont.dsea) * ff_object->fragmentation(DSSpiNLO::down, hadron);
    qqfactor += pdf_object->cont.dsea * ff_object->fragmentation(DSSpiNLO::down_bar, hadron);
    qqfactor += pdf_object->cont.str * ff_object->fragmentation(DSSpiNLO::strange, hadron);
    qqfactor += pdf_object->cont.sbar * ff_object->fragmentation(DSSpiNLO::strange_bar, hadron);
    
    if (ctx->projectile == deuteron) {
        // Neutron contributions (for deuteron collisions), assuming isospin symmetry:
        qqfactor += (pdf_object->cont.dnv + pdf_object->cont.dsea) * ff_object->fragmentation(DSSpiNLO::up, hadron);
        qqfactor += pdf_object->cont.dsea * ff_object->fragmentation(DSSpiNLO::up_bar, hadron);
        qqfactor += (pdf_object->cont.upv + pdf_object->cont.usea) * ff_object->fragmentation(DSSpiNLO::down, hadron);
        qqfactor += pdf_object->cont.usea * ff_object->fragmentation(DSSpiNLO::down_bar, hadron);
        qqfactor += pdf_object->cont.str * ff_object->fragmentation(DSSpiNLO::strange, hadron);
        qqfactor += pdf_object->cont.sbar * ff_object->fragmentation(DSSpiNLO::strange_bar, hadron);
    }
    
    this->qqfactor = qqfactor;

    
    // Proton contribution:
    ggfactor += pdf_object->cont.glu * ff_object->fragmentation(DSSpiNLO::gluon, hadron);
    
    if (ctx->projectile == deuteron) {
        // Neutron contribution (for deuteron collisions), assuming isospin symmetry:
        ggfactor *= 2;
    }
    
    this->ggfactor = ggfactor;

    
    // Proton contributions:
    gqfactor += (  pdf_object->cont.upv + 2*pdf_object->cont.usea
                 + pdf_object->cont.dnv + 2*pdf_object->cont.dsea
                 + pdf_object->cont.str
                 + pdf_object->cont.sbar
                ) * ff_object->fragmentation(DSSpiNLO::gluon, hadron);
    
    if (ctx->projectile == deuteron) {
        // Neutron contributions (for deuteron collisions), assuming isospin symmetry:
        gqfactor *= 2;
    }
    
    this->gqfactor = gqfactor;

    
    // Proton contributions:
    qgfactor += pdf_object->cont.glu * (  ff_object->fragmentation(DSSpiNLO::up, hadron)
                                        + ff_object->fragmentation(DSSpiNLO::up_bar, hadron)
                                        + ff_object->fragmentation(DSSpiNLO::down, hadron)
                                        + ff_object->fragmentation(DSSpiNLO::down_bar, hadron)
                                        + ff_object->fragmentation(DSSpiNLO::strange, hadron)
                                        + ff_object->fragmentation(DSSpiNLO::strange_bar, hadron));
    
    if (ctx->projectile == deuteron) {
        // Neutron contributions (for deuteron collisions), assuming isospin symmetry:
        qgfactor *= 2;
    }
    
    this->qgfactor = qgfactor;
}

void IntegrationContext::update_positions(double xx, double xy, double yx, double yy, double bx, double by) {
    checkfinite(xx);
    checkfinite(xy);
    checkfinite(yx);
    checkfinite(yy);
    checkfinite(bx);
    checkfinite(by);
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

    // Calculate the new gluon distribution values
    // this has to be done after kinematics are updated
    if (r2 > 0) {
        this->S2r = ctx->gdist->S2(r2, this->Yg);
        if (s2 > 0 || t2 > 0) {
            this->S4rst = ctx->gdist->S4(r2, s2, t2, this->Yg);
        }
        else {
            this->S4rst = NAN;
        }
    }
    else {
        this->S2r = NAN;
    }
}

void IntegrationContext::update_momenta(double q1x, double q1y, double q2x, double q2y, double q3x, double q3y) {
    checkfinite(q1x);
    checkfinite(q1y);
    checkfinite(q2x);
    checkfinite(q2y);
    checkfinite(q3x);
    checkfinite(q3y);
    this->q1x = q1x;
    this->q1y = q1y;
    this->q2x = q2x;
    this->q2y = q2y;
    this->q3x = q3x;
    this->q3y = q3y;
    this->q12 = q1x*q1x + q1y*q1y;
    this->q22 = q2x*q2x + q2y*q2y;
    this->q32 = q3x*q3x + q3y*q3y;
    
    if (q12 > 0) {
        this->Fq1 = ctx->gdist->F(q12, this->Yg);
        this->Fkq1 = ctx->gdist->F((kT - q1x)*(kT - q1x) + q1y*q1y, this->Yg);
    }
    else {
        this->Fq1 = NAN;
        this->Fkq1 = NAN;
    }
    if (q22 > 0) {
        this->Fq2 = ctx->gdist->F(q22, this->Yg);
        this->Fkq2 = ctx->gdist->F((kT - q2x)*(kT - q2x) + q2y*q2y, this->Yg);
    }
    else {
        this->Fq2 = NAN;
        this->Fkq2 = NAN;
    }
    if (q32 > 0) {
        this->Fq3 = ctx->gdist->F(q32, this->Yg);
        this->Fkq3 = ctx->gdist->F((kT - q3x)*(kT - q3x) + q3y*q3y, this->Yg);
    }
    else {
        this->Fq3 = NAN;
        this->Fkq3 = NAN;
    }
}

void IntegrationContext::update_auxiliary(double xiprime) {
    checkfinite(xiprime);
    assert(xiprime >= 0 && xiprime <= 1);
    this->xiprime = xiprime;
    this->xiprime2 = xiprime * xiprime;
}

