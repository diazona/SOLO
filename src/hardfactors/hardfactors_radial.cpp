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
#include <gsl/gsl_sf_bessel.h>
#include "hardfactor.h"
#include "hardfactors_radial.h"

#define checkfinite(d) assert(gsl_finite(d))

namespace radial {
    const AngleIndependentPositionIntegrationType dipole(1);
    const AngleIndependentPositionIntegrationType quadrupole(2);
    const RescaledAngleIndependentPositionIntegrationType rescaled_dipole(1);
    const RescaledAngleIndependentPositionIntegrationType rescaled_quadrupole(2);
}

using namespace radial;

void H02qq::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = 0.25*M_1_PI*M_1_PI * ictx->ctx->Sperp / ictx->z2 * ictx->qqfactor * ictx->S2r;
    double r = sqrt(ictx->r2);
    double b_phase = ictx->kT * r;
    *real = amplitude * gsl_sf_bessel_J0(b_phase);
    *imag = 0;
}

void H12qq::Fs(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx->ctx->c0r_optimization) {
        *real = *imag = 0;
        return;
    }
    double amplitude = 0.25*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->CF * ictx->ctx->Sperp / ictx->z2 * ictx->qqfactor * ictx->S2r *
        (1 + ictx->xi2) * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
    double r = sqrt(ictx->r2);
    double b_phase1 = ictx->kT * r;
    double b_phase2 = ictx->kT * r / ictx->xi;
    *real = amplitude * (gsl_sf_bessel_J0(b_phase1) + gsl_sf_bessel_J0(b_phase2) / ictx->xi2);
    *imag = 0;
}

void H12qq::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = 0.75*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->CF * ictx->ctx->Sperp / ictx->z2 * ictx->qqfactor * ictx->S2r;
    assert(ictx->xi == 1);
    assert(ictx->xi2 == 1);
    if (ictx->ctx->c0r_optimization) {
        double term2 = -(-2*M_EULER + log(4) - log(ictx->r2 * ictx->kT2));
        double r = sqrt(ictx->r2);
        double phase1 = ictx->kT * r;
        *real = amplitude * term2 * gsl_sf_bessel_J0(phase1);
        *imag = 0;
    }
    else {
        // this expression can be simplified by assuming xi=1 and combining
        // terms, but I don't do that so that it is obvious how the expressions
        // here match up to the formulas in the paper
        double term1 = 0.5 * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
        double term2 = -(-2*M_EULER + log(4) - log(ictx->r2 * ictx->kT2));
        double r = sqrt(ictx->r2);
        double phase1 = ictx->kT * r;
        double phase2 = ictx->kT * r / ictx->xi;
        *real = amplitude * (term1 * (gsl_sf_bessel_J0(phase1) + gsl_sf_bessel_J0(phase2) / ictx->xi2) + term2 * gsl_sf_bessel_J0(phase1));
        *imag = 0;
    }
}

void H12qq1::Fs(const IntegrationContext* ictx, double* real, double* imag) const {
    // implementation is the same as H12qq:Fs
    if (ictx->ctx->c0r_optimization) {
        *real = *imag = 0;
        return;
    }
    double amplitude = 0.25*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->CF * ictx->ctx->Sperp / ictx->z2 * ictx->qqfactor * ictx->S2r *
        (1 + ictx->xi2) * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
    double r = sqrt(ictx->r2);
    double b_phase1 = ictx->kT * r;
    double b_phase2 = ictx->kT * r / ictx->xi;
    *real = amplitude * (gsl_sf_bessel_J0(b_phase1) + gsl_sf_bessel_J0(b_phase2) / ictx->xi2);
    *imag = 0;
}

void H12qq1A::Fs(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx->ctx->c0r_optimization) {
        *real = *imag = 0;
        return;
    }
    double amplitude = 0.25*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->CF * ictx->ctx->Sperp / ictx->z2 * ictx->qqfactor * ictx->S2r *
        (1 + ictx->xi2) * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
    double r = sqrt(ictx->r2);
    double b_phase1 = ictx->kT * r;
    *real = amplitude * gsl_sf_bessel_J0(b_phase1);
    *imag = 0;
}

void H12qq1B::Fs(const IntegrationContext* ictx, double* real, double* imag) const {
    // implementation is the same as H12qq:Fs
    if (ictx->ctx->c0r_optimization) {
        *real = *imag = 0;
        return;
    }
    double amplitude = 0.25*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->CF * ictx->ctx->Sperp / ictx->z2 * ictx->qqfactor * ictx->S2r *
        (1 + ictx->xi2) * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
    double r = sqrt(ictx->r2);
    double b_phase2 = ictx->kT * r / ictx->xi;
    *real = amplitude * gsl_sf_bessel_J0(b_phase2) / ictx->xi2;
    *imag = 0;
}

void H12qq2::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    // implemented the same as H12qq:Fd setting term2 to zero
    double amplitude = 0.75*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->CF * ictx->ctx->Sperp / ictx->z2 * ictx->qqfactor * ictx->S2r;
    assert(ictx->xi == 1);
    assert(ictx->xi2 == 1);
    if (ictx->ctx->c0r_optimization) {
        *real = 0;
        *imag = 0;
    }
    else {
        double term1 = 0.5 * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
        double r = sqrt(ictx->r2);
        double phase1 = ictx->kT * r;
        double phase2 = ictx->kT * r / ictx->xi;
        *real = amplitude * term1 * (gsl_sf_bessel_J0(phase1) + gsl_sf_bessel_J0(phase2) / ictx->xi2);
        *imag = 0;
    }
}

void H12qq3::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    // implemented the same as H12qq:Fd setting term1 to zero
    double amplitude = 0.75*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->CF * ictx->ctx->Sperp / ictx->z2 * ictx->qqfactor * ictx->S2r;
    assert(ictx->xi == 1);
    assert(ictx->xi2 == 1);
    double term2 = -(-2*M_EULER + log(4) - log(ictx->r2 * ictx->kT2));
    double r = sqrt(ictx->r2);
    double phase1 = ictx->kT * r;
    *real = amplitude * term2 * gsl_sf_bessel_J0(phase1);
    *imag = 0;
}


void H012qqExp::Fs(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx->ctx->c0r_optimization) {
        *real = *imag = 0;
        return;
    }
    double amplitude = 0.25*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->CF * ictx->ctx->Sperp / ictx->z2 * ictx->qqfactor * ictx->S2r *
        (1 + ictx->xi2) * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
    double r = sqrt(ictx->r2);
    double b_phase1 = ictx->kT * r;
    double b_phase2 = ictx->kT * r / ictx->xi;
    *real = amplitude * (gsl_sf_bessel_J0(b_phase1) + gsl_sf_bessel_J0(b_phase2) / ictx->xi2);
    *imag = 0;
}

void H012qqExp::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = 0.25*M_1_PI*M_1_PI * ictx->ctx->Sperp / ictx->z2 * ictx->qqfactor * ictx->S2r;
    double factor1 = 3 * ictx->alphas_2pi * ictx->ctx->CF;
    double r = sqrt(ictx->r2);
    assert(ictx->xi == 1);
    assert(ictx->xi2 == 1);
    if (ictx->ctx->c0r_optimization) {
        double term2 = -(-2*M_EULER + log(4) - log(ictx->r2 * ictx->kT2));
        double b_phase = ictx->kT * r;
        *real = amplitude * exp(factor1 * term2);
        *imag = 0;
    }
    else {
        double logterm = log(ictx->mu2 / ictx->kT2);
        double b_phase = ictx->kT * r;
        *real = amplitude * exp(factor1 * logterm) * gsl_sf_bessel_J0(b_phase);
        *imag = 0;
    }
}

void H12gg::Fs(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx->ctx->c0r_optimization) {
        *real = *imag = 0;
        return;
    }
    double amplitude = 1/(4*M_PI*M_PI) * ictx->alphas_2pi * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r *
        ictx->ctx->Nc * 2 * ictx->xi
        * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
    double r = sqrt(ictx->r2);
    double b_phase1 = ictx->kT * r;
    double b_phase2 = ictx->kT * r / ictx->xi;
    *real = amplitude * (gsl_sf_bessel_J0(b_phase1) + gsl_sf_bessel_J0(b_phase2) / ictx->xi2);
    *imag = 0;
}

void H12gg::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx->ctx->c0r_optimization) {
        *real = *imag = 0;
        return;
    }
    double amplitude = 1/(4*M_PI*M_PI) * ictx->alphas_2pi * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r *
        ictx->ctx->Nc * 2 * (1.0 / ictx->xi - 1.0 + ictx->xi - ictx->xi2)
        * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
    double r = sqrt(ictx->r2);
    double b_phase1 = ictx->kT * r;
    double b_phase2 = ictx->kT * r / ictx->xi;
    *real = amplitude * (gsl_sf_bessel_J0(b_phase1) + gsl_sf_bessel_J0(b_phase2) / ictx->xi2);
    *imag = 0;
}

void H12gg::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx->ctx->c0r_optimization) {
        double amplitude = 1/(4*M_PI*M_PI) * ictx->alphas_2pi * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r;
        double term2 = -ictx->ctx->Nc * 2 * (11.0/6.0 - 2 * ictx->ctx->Nf * ictx->ctx->TR / (3 * ictx->ctx->Nc)) // independent term
            * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->kT2));
        double r = sqrt(ictx->r2);
        double b_phase1 = ictx->kT * r;
        *real = amplitude * term2 * gsl_sf_bessel_J0(b_phase1);
        *imag = 0;
    }
    else {
        double amplitude = 1/(4*M_PI*M_PI) * ictx->alphas_2pi * ictx->ggfactor / ictx->z2 * ictx->S2r * ictx->S2r;
        double term1 = ictx->ctx->Nc * (11.0/6.0 - 2 * ictx->ctx->Nf * ictx->ctx->TR / (3 * ictx->ctx->Nc))       // P_gg term
            * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
        double term2 = -ictx->ctx->Nc * 2 * (11.0/6.0 - 2 * ictx->ctx->Nf * ictx->ctx->TR / (3 * ictx->ctx->Nc)) // independent term
            * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->kT2));
        double r = sqrt(ictx->r2);
        double b_phase1 = ictx->kT * r;
        double b_phase2 = ictx->kT * r / ictx->xi;
        double phase1 = ictx->kT * ictx->rx;
        double phase2 = ictx->kT * ictx->rx / ictx->xi;
        *real = amplitude * (term1 * (gsl_sf_bessel_J0(b_phase1) + gsl_sf_bessel_J0(b_phase2) / ictx->xi2) + term2 * gsl_sf_bessel_J0(b_phase1));
        *imag = 0;
    }
}

void H112gq::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx->ctx->c0r_optimization) {
        *real = *imag = 0;
        return;
    }
    double amplitude = 1/(4*M_PI*M_PI) * ictx->alphas_2pi * ictx->gqfactor / ictx->z2 * ictx->S2r *
        0.5 * ictx->ctx->Nc * (1 + (1 - ictx->xi)*(1 - ictx->xi)) / (ictx->xi * ictx->xi2) * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
    double r = sqrt(ictx->r2);
    double b_phase = ictx->kT * r / ictx->xi;
    checkfinite(amplitude);
    checkfinite(b_phase);
    *real = amplitude * gsl_sf_bessel_J0(b_phase);
    *imag = 0;
}

void H122gq::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    if (ictx->ctx->c0r_optimization) {
        *real = *imag = 0;
        return;
    }
    double amplitude = 1/(4*M_PI*M_PI) * ictx->alphas_2pi * ictx->gqfactor / ictx->z2 * ictx->S2r * ictx->S2r *
        0.5 * ictx->ctx->Nc * (1 + (1 - ictx->xi)*(1 - ictx->xi)) / ictx->xi * (-2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2));
    double r = sqrt(ictx->r2);
    double b_phase = ictx->kT * r;
    checkfinite(amplitude);
    checkfinite(b_phase);
    *real = amplitude * gsl_sf_bessel_J0(b_phase);
    *imag = 0;
}

void H112qg::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = 0.125*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->Sperp / ictx->z2 * ictx->qgfactor * ictx->S2r *
        (1 - 2 * ictx->xi + 2 * ictx->xi2) * (ictx->ctx->c0r_optimization ? -1 : -2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2) - 1);
    double r = sqrt(ictx->r2);
    double b_phase = ictx->kT * r;
    checkfinite(amplitude);
    checkfinite(b_phase);
    *real = amplitude * gsl_sf_bessel_J0(b_phase);
    *imag = 0;
}

void H122qg::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    double amplitude = 0.125*M_1_PI*M_1_PI * ictx->alphas_2pi * ictx->ctx->Sperp / ictx->z2 * ictx->qgfactor * ictx->S2r * ictx->S2r *
        (1 / ictx->xi2 - 2 / ictx->xi + 2) * (ictx->ctx->c0r_optimization ? -1 : -2*M_EULER + log(4) - log(ictx->r2 * ictx->mu2) - 1);
    double r = sqrt(ictx->r2);
    double b_phase = ictx->kT * r / ictx->xi;
    checkfinite(amplitude);
    checkfinite(b_phase);
    *real = amplitude * gsl_sf_bessel_J0(b_phase);
    *imag = 0;
}

