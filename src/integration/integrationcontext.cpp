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

void IntegrationContext::recalculate_from_position(const bool quadrupole) {
    checkfinite(xx);
    checkfinite(xy);
    checkfinite(yx);
    checkfinite(yy);
    rx = xx - yx;
    ry = xy - yy;
    r2 = rx*rx + ry*ry;
    if (quadrupole) {
        checkfinite(bx);
        checkfinite(by);
        sx = xx - bx;
        sy = xy - by;
        s2 = sx*sx + sy*sy;
        tx = yx - bx;
        ty = yy - by;
        t2 = tx*tx + ty*ty;
    }
}

void IntegrationContext::recalculate_from_momentum(const size_t dimensions) {
    switch (dimensions) { // intentionally omitting break statements
        case 3:
            checkfinite(q3x);
            checkfinite(q3y);
            q32 = q3x*q3x + q3y*q3y;
            assert(q32 >= 0);
        case 2:
            checkfinite(q2x);
            checkfinite(q2y);
            q22 = q2x*q2x + q2y*q2y;
            assert(q22 >= 0);
        case 1:
            checkfinite(q1x);
            checkfinite(q1y);
            q12 = q1x*q1x + q1y*q1y;
            assert(q12 >= 0);
        case 0:
            break;
        default:
            assert(false);
    }
}

void IntegrationContext::recalculate_longitudinal() {
    z2 = z*z;
    kT2 = ctx.pT2 / z2;
    kT = sqrt(kT2);
    xi2 = xi * xi;
    /* NOTE: this conflicts with the xp used by previous versions of SOLO.
     * Previous versions defined it as tau / (z * xi). We now express this
     * as xp / xi, consistent with the notation in e.g. arxiv:1604.00225
     */
    xp = ctx.tau / z;
    xg = kT / ctx.sqs * exp(-ctx.Y);
    Yg = -log(xg);
}

void IntegrationContext::recalculate_exact_longitudinal() {
    xa = exp(-ctx.Y) / ctx.sqs * (kT + ((kT - q1x) * (kT - q1x) + q1y * q1y) / kT * xi / (1 - xi));
    Ya = -log(xa);
}

void IntegrationContext::recalculate_coupling() {
    alphas = ctx.cpl->alphas(kT2);
    alphas_2pi = alphas * 0.5 * M_1_PI;
}

void IntegrationContext::recalculate_position_gdist(const bool quadrupole) {
    S2r = ctx.gdist->S2(r2, Yg);
    S4rst = quadrupole ? ctx.gdist->S4(r2, s2, t2, Yg) : NAN;
}

void IntegrationContext::recalculate_momentum_gdist(const size_t dimensions, const bool exact) {
    double Y = exact ? Ya : Yg;
    Qs2 = ctx.gdist->Qs2(Y);
    switch (dimensions) { // intentionally omitting break statements
        case 3:
            Fq3 = ctx.gdist->F(q32, Y);
            Fkq3 = ctx.gdist->F((kT - q3x) * (kT - q3x) + q3y * q3y, Y);
        case 2:
            Fq2 = ctx.gdist->F(q22, Y);
            Fkq2 = ctx.gdist->F((kT - q2x) * (kT - q2x) + q2y * q2y, Y);
        case 1:
            Fq1 = ctx.gdist->F(q12, Y);
            Fkq1 = ctx.gdist->F((kT - q1x) * (kT - q1x) + q1y * q1y, Y);
        case 0:
            Fk = ctx.gdist->F(kT2, Y);
            break;
        default:
            assert(false);
    }
}

void IntegrationContext::recalculate_parton_functions(const bool divide_xi) {
    // Finally, update the parton functions
    double qqfactor = 0.0, ggfactor = 0.0, gqfactor = 0.0, qgfactor = 0.0;

    c_mstwpdf* pdf_object = tlctx.pdf_object;
    DSSpiNLO* ff_object = tlctx.ff_object;

    DSSpiNLO::hadron hadron = ctx.hadron;

    // Calculate the new quark/gluon factors
    mu2 = ctx.fs->mu2(*this);
    pdf_object->update(divide_xi ? xp / xi : xp, sqrt(mu2));
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

void IntegrationContext::recalculate_everything_from_position(const bool quadrupole, const bool divide_xi) {
    recalculate_from_position(quadrupole);
    recalculate_longitudinal();
    recalculate_coupling();
    recalculate_position_gdist(quadrupole);
    recalculate_parton_functions(divide_xi);
}

void IntegrationContext::recalculate_everything_from_momentum(const size_t dimensions, const bool exact_xg, const bool divide_xi) {
    recalculate_from_momentum(dimensions);
    recalculate_longitudinal();
    if (exact_xg) {
        recalculate_exact_longitudinal();
    }
    recalculate_coupling();
    recalculate_momentum_gdist(dimensions, exact_xg);
    recalculate_parton_functions(divide_xi);
}

void IntegrationContext::recalculate_everything(const bool exact_xg, const bool divide_xi) {
    recalculate_from_position(true);
    recalculate_from_momentum(3);
    recalculate_longitudinal();
    recalculate_exact_longitudinal();
    recalculate_coupling();
    // use some heuristics to guess whether to compute position or momentum gluon distributions
    if (r2 > 0 || s2 > 0 || t2 > 0) {
        recalculate_position_gdist(s2 > 0 || t2 > 0);
    }
    if (q12 > 0 || q22 > 0 || q32 > 0) {
        recalculate_momentum_gdist(3, exact_xg);
    }
    else {
        recalculate_momentum_gdist(0);
    }
    recalculate_parton_functions(divide_xi);
}
