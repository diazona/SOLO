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

// Hard factors in position space

#include <cassert>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include "hardfactor.h"
#include "hardfactors_position.h"

#define checkfinite(d) assert(gsl_finite(d))

using namespace position;

void H02qq::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = 0.25*M_1_PI*M_1_PI * ictx->ctx->Sperp / ictx->z2 * ictx->qqfactor * ictx->S2r;
    double phase = -ictx->kT * ictx->rx; // take angle of k_perp to be 0
    *real = amplitude * cos(phase);
    *imag = amplitude * sin(phase);
}

void H12qq::Fs(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx->ctx->c0r_optimization) {
        *real = *imag = 0;
        return;
    }
    double amplitude = 0.25*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->CF * ictx->ctx->Sperp / ictx->z2 * ictx->qqfactor * ictx->S2r *
        (1 + ictx->xi2) * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
    double phase1 = -ictx->kT * ictx->rx;
    double phase2 = -ictx->kT * ictx->rx / ictx->xi;
    *real = amplitude * (cos(phase1) + cos(phase2) / ictx->xi2);
    *imag = amplitude * (sin(phase1) + sin(phase2) / ictx->xi2);
}
void H12qq::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = 0.75*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->CF * ictx->ctx->Sperp / ictx->z2 * ictx->qqfactor * ictx->S2r;
    if (ictx->ctx->c0r_optimization) {
        double term2 = -(-2*M_EULER + log(4) - log(ictx->r2 * ictx->kT2));
        double phase1 = -ictx->kT * ictx->rx;
        *real = amplitude * term2 * cos(phase1);
        *imag = amplitude * term2 * sin(phase1);
    }
    else {
        double term1 = 0.5 * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
        double term2 = -(-2*M_EULER + log(4) - log(ictx->r2 * ictx->kT2));
        double phase1 = -ictx->kT * ictx->rx;
        double phase2 = -ictx->kT * ictx->rx / ictx->xi;
        *real = amplitude * (term1 * (cos(phase1) + cos(phase2) / ictx->xi2) + term2 * cos(phase1));
        *imag = amplitude * (term1 * (sin(phase1) + sin(phase2) / ictx->xi2) + term2 * sin(phase1));
    }
}

void I1(double x, double* real, double* imag) {
    bool negative = false;
    if (x < 0) {
        // if argument is negative, switch it to positive and flip the sign of the imaginary part of the result
        x = -x;
        negative = true;
    }
    *real = -2 * M_EULER - gsl_sf_sinc(x * M_1_PI) + (cos(x) - 1) / (x * x) + 2 * gsl_sf_Ci(x) - 2 * log(x);
    *imag = (negative ? -1 : 1) * ((2 - cos(x) - gsl_sf_sinc(x * M_1_PI)) / x - 2 * gsl_sf_Si(x));
}

void H14qqPrimary::Fs(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = -0.25*M_1_PI*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->Nc * ictx->ctx->Sperp / ictx->z2 * ictx->qqfactor * ictx->S4rst
        * (1 + ictx->xi2) / ictx->xi
        * (ictx->sx * ictx->tx + ictx->sy * ictx->ty)
        / ( ictx->s2 * ictx->t2 );
    double phase = -ictx->kT * (ictx->sx / ictx->xi - ictx->tx);
    *real = amplitude * cos(phase);
    *imag = amplitude * sin(phase);
}

void H14qqPrimary::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double realI, imagI;
    I1(ictx->kT * ictx->tx, &realI, &imagI);
    double amplitude = 0.25*M_1_PI*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->Nc * ictx->ctx->Sperp / ictx->z2 * ictx->qqfactor
        * (ictx->S4rst - ictx->ctx->gdist->S4(ictx->s2, ictx->s2, 0, ictx->Yg))
        / ictx->t2;
    double phase = -ictx->kT * ictx->rx;
    *real = amplitude * (cos(phase) * realI - sin(phase) * imagI);
    *imag = amplitude * (cos(phase) * imagI + sin(phase) * realI);
}

void H14qqResidual::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = (1.0/6.0 - 0.625*M_1_PI*M_1_PI) * ictx->alphas_2pi * ictx->ctx->Nc * ictx->ctx->Sperp / ictx->z2 * ictx->qqfactor
        * ictx->ctx->gdist->S4(ictx->r2, ictx->r2, 0, ictx->Yg);
    double phase = ictx->kT * ictx->rx;
    *real = amplitude * cos(phase);
    *imag = amplitude * sin(phase);
}

void H02gg::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = 1/(4*M_PI*M_PI) * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r;
    double phase = -ictx->kT * ictx->rx;
    *real = amplitude * cos(phase);
    *imag = amplitude * sin(phase);
}

void H12gg::Fs(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx->ctx->c0r_optimization) {
        *real = *imag = 0;
        return;
    }
    double amplitude = 1/(4*M_PI*M_PI) * ictx->alphas_2pi * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r *
        ictx->ctx->Nc * 2 * ictx->xi
        * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
    double phase1 = -ictx->kT * ictx->rx;
    double phase2 = -ictx->kT * ictx->rx / ictx->xi;
    *real = amplitude * (cos(phase1) + cos(phase2) / ictx->xi2);
    *imag = amplitude * (sin(phase1) + sin(phase2) / ictx->xi2);
}

void H12gg::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx->ctx->c0r_optimization) {
        *real = *imag = 0;
        return;
    }
    double amplitude = 1/(4*M_PI*M_PI) * ictx->alphas_2pi * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r *
        ictx->ctx->Nc * 2 * (1.0 / ictx->xi - 1.0 + ictx->xi - ictx->xi2)
        * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
    double phase1 = ictx->kT * ictx->rx;
    double phase2 = ictx->kT * ictx->rx / ictx->xi;
    *real = amplitude * (cos(phase1) + cos(phase2) / ictx->xi2);
    *imag = amplitude * (sin(phase1) + sin(phase2) / ictx->xi2);
}

void H12gg::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = 1/(4*M_PI*M_PI) * ictx->alphas_2pi * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r;
    if (ictx->ctx->c0r_optimization) {
        double term2 = -ictx->ctx->Nc * 2 * (11.0/6.0 - 2 * ictx->ctx->Nf * ictx->ctx->TR / (3 * ictx->ctx->Nc)) // independent term
            * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->kT2));
        double phase1 = ictx->kT * ictx->rx;
        *real = amplitude * term2 * cos(phase1);
        *imag = amplitude * term2 * cos(phase1);
    }
    else {
        double term1 = ictx->ctx->Nc * (11.0/6.0 - 2 * ictx->ctx->Nf * ictx->ctx->TR / (3 * ictx->ctx->Nc))       // P_gg term
            * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
        double term2 = -ictx->ctx->Nc * 2 * (11.0/6.0 - 2 * ictx->ctx->Nf * ictx->ctx->TR / (3 * ictx->ctx->Nc)) // independent term
            * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->kT2));
        double phase1 = ictx->kT * ictx->rx;
        double phase2 = ictx->kT * ictx->rx / ictx->xi;
        *real = amplitude * (term1 * (cos(phase1) + cos(phase2) / ictx->xi2) + term2 * cos(phase1));
        *imag = amplitude * (term1 * (sin(phase1) + sin(phase2) / ictx->xi2) + term2 * cos(phase1));
    }
}

void I2(double x, double* real, double* imag) {
    double sx = sin(x), cx = cos(x);
    double ix = 1.0 / x;
    double ix2 = ix * ix;
    double ix3 = ix2 * ix;
    *real = -4 * ix3 * sx + 2 * ix2 * (1 + cx) + ix * sx;
    *imag = 4 * ix3 * (1 - cx) - 2 * ix2 * sx - ix * (1 - cx);
}

void H12qqbarPrimary::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double realI, imagI;
    I2(ictx->kT * sqrt(ictx->r2), &realI, &imagI);
    double amplitude = 1/(4*M_PI*M_PI*4*M_PI*M_PI) * ictx->alphas_2pi * ictx->ggfactor / ictx->z2 *
        8 * M_PI * ictx->ctx->Nf * ictx->ctx->TR * ictx->ctx->Sperp
        * (ictx->ctx->gdist->S2(ictx->s2, ictx->Yg) - ictx->ctx->gdist->S2(ictx->t2, ictx->Yg)) * ictx->ctx->gdist->S2(ictx->t2, ictx->Yg)
        / ictx->r2;
    double phase = ictx->kT * ictx->tx;
    *real = amplitude * (cos(phase) * realI + sin(phase) * imagI);
    *imag = amplitude * (cos(phase) * imagI - sin(phase) * realI);
}


void H12qqbarResidual::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = 1/(4*M_PI*M_PI) * ictx->alphas_2pi * ictx->ggfactor / ictx->z2 * ictx->ctx->Sperp * ictx->S2r * ictx->S2r *
        4 * M_PI * ictx->ctx->Nf * ictx->ctx->TR * (13.0 / 18.0);
    double phase = -ictx->kT * ictx->tx;
    *real = amplitude * cos(phase);
    *imag = amplitude * sin(phase);
}

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

void H16ggPrimary::Fs(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = M_1_PI*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->Nc * ictx->ctx->Sperp / ictx->z2
        * ictx->ggfactor * ictx->ctx->gdist->S2(ictx->s2, ictx->Yg) * ictx->ctx->gdist->S2(ictx->t2, ictx->Yg) * ictx->S2r * (1 - ictx->xi + ictx->xi2)*(1 - ictx->xi + ictx->xi2) / ictx->xi2
        * (ictx->rx*ictx->tx + ictx->ry*ictx->ty) / (ictx->r2 * ictx->t2);
    double phase = -ictx->kT * (ictx->rx + ictx->tx / ictx->xi);
    *real = amplitude * cos(phase);
    *imag = amplitude * sin(phase);
}

void H16ggPrimary::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double realI, imagI;
    I3(ictx->kT * ictx->tx, &realI, &imagI);
    double amplitude = M_1_PI*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->Nc * ictx->ctx->Sperp / ictx->z2
        * ictx->ggfactor * ictx->ctx->gdist->S2(ictx->s2, ictx->Yg) * (ictx->ctx->gdist->S2(ictx->t2, ictx->Yg) * ictx->ctx->gdist->S2((ictx->sx + ictx->tx)*(ictx->sx + ictx->tx) + (ictx->sy + ictx->ty)*(ictx->sy + ictx->ty), ictx->Yg) - ictx->ctx->gdist->S2(ictx->s2, ictx->Yg) * ictx->ctx->gdist->S2(0, ictx->Yg))
        / ictx->t2;
    double phase = -ictx->kT * ictx->rx;
    *real = amplitude * (cos(phase) * realI - sin(phase) * imagI);
    *imag = amplitude * (cos(phase) * imagI + sin(phase) * realI);
}

void H16ggResidual::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = (-67.0 / 36.0 * M_1_PI*M_1_PI + 1.0 / 3.0) * ictx->alphas_2pi * ictx->ctx->Nc * ictx->ctx->Sperp / ictx->z2
        * ictx->ggfactor * ictx->S2r * ictx->S2r * ictx->ctx->gdist->S2(0, ictx->Yg);
    double phase = -ictx->kT * ictx->rx;
    *real = amplitude * cos(phase);
    *imag = amplitude * cos(phase);
}

void H112gq::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx->ctx->c0r_optimization) {
        *real = *imag = 0;
        return;
    }
    double amplitude = 1/(4*M_PI*M_PI) * ictx->alphas_2pi * ictx->gqfactor / ictx->z2 * ictx->S2r *
        0.5 * ictx->ctx->Nc * (1 + (1 - ictx->xi)*(1 - ictx->xi)) / (ictx->xi * ictx->xi2) * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
    double phase = -ictx->kT * ictx->rx / ictx->xi;
    checkfinite(amplitude);
    checkfinite(phase);
    *real = amplitude * cos(phase);
    *imag = amplitude * sin(phase);
}

void H122gq::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx->ctx->c0r_optimization) {
        *real = *imag = 0;
        return;
    }
    double amplitude = 1/(4*M_PI*M_PI) * ictx->alphas_2pi * ictx->gqfactor / ictx->z2 * ictx->S2r * ictx->S2r *
        0.5 * ictx->ctx->Nc * (1 + (1 - ictx->xi)*(1 - ictx->xi)) / ictx->xi * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
    double phase = -ictx->kT * ictx->rx;
    checkfinite(amplitude);
    checkfinite(phase);
    *real = amplitude * cos(phase);
    *imag = amplitude * sin(phase);
}

void H14gq::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = 1/(4*M_PI*M_PI*4*M_PI*M_PI) * ictx->alphas_2pi * ictx->gqfactor / ictx->z2 * ictx->S4rst *
        4 * M_PI * ictx->ctx->Nc * (1 + (1 - ictx->xi)*(1 - ictx->xi)) / (ictx->xi2)
        * (ictx->rx*ictx->tx + ictx->ry*ictx->ty)
        / ( ictx->r2 * ictx->t2 );
    double phase = -(ictx->kT * ictx->rx / ictx->xi + ictx->kT * ictx->tx);
    checkfinite(amplitude);
    checkfinite(phase);
    *real = amplitude * cos(phase);
    *imag = amplitude * sin(phase);
}

void H112qg::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = 0.125*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->Sperp / ictx->z2 * ictx->qgfactor * ictx->S2r *
        (1 - 2 * ictx->xi + 2 * ictx->xi2) * (ictx->ctx->c0r_optimization ? -1 : -2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2) - 1);
    double phase = -ictx->kT * ictx->rx;
    checkfinite(amplitude);
    checkfinite(phase);
    *real = amplitude * cos(phase);
    *imag = amplitude * sin(phase);
}

void H122qg::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = 0.125*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->Sperp / ictx->z2 * ictx->qgfactor * ictx->S2r * ictx->S2r *
        (1 / ictx->xi2 - 2 / ictx->xi + 2) * (ictx->ctx->c0r_optimization ? -1 : -2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2) - 1);
    double phase = -ictx->kT * ictx->rx / ictx->xi;
    checkfinite(amplitude);
    checkfinite(phase);
    *real = amplitude * cos(phase);
    *imag = amplitude * sin(phase);
}

void H14qg::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = 0.25*M_1_PI*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->Sperp / ictx->z2 * ictx->qgfactor * ictx->S4rst *
        (1 / ictx->xi - 2 + 2 * ictx->xi)
        * (ictx->rx*ictx->tx + ictx->ry*ictx->ty)
        / ( ictx->r2 * ictx->t2 );
    double phase = -(ictx->kT * ictx->rx + ictx->kT * ictx->tx / ictx->xi);
    checkfinite(amplitude);
    checkfinite(phase);
    *real = amplitude * cos(phase);
    *imag = amplitude * sin(phase);
}

