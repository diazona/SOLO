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
#include "../mstwpdf.h"
#include "../dss_pinlo/dss_pinlo.h"
#include "integrationcontext.h"
#include "../gluondist/gluondist.h"

#define checkfinite(d) assert(gsl_finite(d))

void IntegrationContext::recalculate_everything() {
    // These are the things that should have been manually set
    checkfinite(z);
    checkfinite(xi);
    checkfinite(xx);
    checkfinite(xy);
    checkfinite(yx);
    checkfinite(yy);
    checkfinite(bx);
    checkfinite(by);
    checkfinite(q1x);
    checkfinite(q1y);
    checkfinite(q2x);
    checkfinite(q2y);
    checkfinite(q3x);
    checkfinite(q3y);

    // Update dipole sizes from positions
    rx = xx - yx;
    ry = xy - yy;
    r2 = rx*rx + ry*ry;
    sx = xx - bx;
    sy = xy - by;
    s2 = sx*sx + sy*sy;
    tx = yx - bx;
    ty = yy - by;
    t2 = tx*tx + ty*ty;

    // Update momenta
    q12 = q1x*q1x + q1y*q1y;
    q22 = q2x*q2x + q2y*q2y;
    q32 = q3x*q3x + q3y*q3y;
    assert(q12 >= 0);
    assert(q22 >= 0);
    assert(q32 >= 0);

    // Update kinematical variables
    z2 = z*z;
    // the definition of kT is also referenced in ExactKinematicIntegrationRegion
    kT2 = ctx.pT2 / z2;
    kT = sqrt(kT2);
    xg = kT / ctx.sqs * exp(-ctx.Y); // doubles as \hat{x}_a
    xi2 = xi*xi;

    // the definition of xp is also referenced in ExactKinematicIntegrationRegion
    xp = ctx.tau / (z * xi);
    xa = exp(-ctx.Y) / ctx.sqs * (kT + ((kT - q1x)*(kT - q1x) + q1y*q1y) / kT * xi / (1 - xi));
    Yg = -log(xg);
    Ya = -log(xa);
    Qs2 = ctx.gdist->Qs2(Ya);
    Fk = ctx.gdist->F(kT2, Ya);
    alphas = ctx.cpl->alphas(kT2);
    alphas_2pi = alphas * 0.5 * M_1_PI;
    xiprime2 = xiprime * xiprime;


    // Calculate the new gluon distribution values
    // this has to be done after kinematics are updated
    if (r2 > 0) {
        S2r = ctx.gdist->S2(r2, Ya);
        if (s2 > 0 || t2 > 0) {
            S4rst = ctx.gdist->S4(r2, s2, t2, Ya);
        }
        else {
            S4rst = NAN;
        }
    }
    else {
        S2r = NAN;
    }

    Fq1 = ctx.gdist->F(q12, Ya);
    Fkq1 = ctx.gdist->F((kT - q1x)*(kT - q1x) + q1y*q1y, Ya);
    Fq2 = ctx.gdist->F(q22, Ya);
    Fkq2 = ctx.gdist->F((kT - q2x)*(kT - q2x) + q2y*q2y, Ya);
    Fq3 = ctx.gdist->F(q32, Ya);
    Fkq3 = ctx.gdist->F((kT - q3x)*(kT - q3x) + q3y*q3y, Ya);

    // Finally, update the parton functions
    double qqfactor = 0.0, ggfactor = 0.0, gqfactor = 0.0, qgfactor = 0.0;

    c_mstwpdf* pdf_object = tlctx.pdf_object;
    DSSpiNLO* ff_object = tlctx.ff_object;

    DSSpiNLO::hadron hadron = ctx.hadron;

    // Calculate the new quark/gluon factors
    mu2 = ctx.fs->mu2(*this);
    pdf_object->update(xp, sqrt(mu2));
    ff_object->update(z, mu2);

    // Proton contributions:
    qqfactor += (pdf_object->cont.upv + pdf_object->cont.usea) * ff_object->fragmentation(DSSpiNLO::up, hadron);
    qqfactor += pdf_object->cont.usea * ff_object->fragmentation(DSSpiNLO::up_bar, hadron);
    qqfactor += (pdf_object->cont.dnv + pdf_object->cont.dsea) * ff_object->fragmentation(DSSpiNLO::down, hadron);
    qqfactor += pdf_object->cont.dsea * ff_object->fragmentation(DSSpiNLO::down_bar, hadron);
    qqfactor += pdf_object->cont.str * ff_object->fragmentation(DSSpiNLO::strange, hadron);
    qqfactor += pdf_object->cont.sbar * ff_object->fragmentation(DSSpiNLO::strange_bar, hadron);

    if (ctx.projectile == deuteron) {
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

    if (ctx.projectile == deuteron) {
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

    if (ctx.projectile == deuteron) {
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

    if (ctx.projectile == deuteron) {
        // Neutron contributions (for deuteron collisions), assuming isospin symmetry:
        qgfactor *= 2;
    }

    this->qgfactor = qgfactor;
}
