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
#include "hardfactors_momentum.h"

#define checkfinite(d) assert(gsl_finite(d))

namespace momentum {
    const NoIntegrationType none;
    const MomentumIntegrationType momentum1(2);
    const MomentumIntegrationType momentum2(4);
    const MomentumIntegrationType momentum3(6);
    const RadialMomentumIntegrationType radialmomentum1(2);
    const RadialMomentumIntegrationType radialmomentum2(4);
    const RadialMomentumIntegrationType radialmomentum3(6);
    const XiPIntegrationType momentumxip1(3);
    const XiPIntegrationType momentumxip2(5);
    const QLimitedMomentumIntegrationType qlim(2);
}

using namespace momentum;

void H02qq::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double value = ictx->ctx->Sperp * ictx->qqfactor / ictx->z2 * ictx->Fk;
    checkfinite(value);
    *real = value;
    *imag = 0;
}

void H14qqSingular::Fs(const IntegrationContext* ictx, double* real, double* imag) const {
    double kA2 = gsl_pow_2(ictx->kT/ictx->xi + ictx->q1x) + gsl_pow_2(ictx->q1y); // (k / xi + q1)^2
    double Fg1 = ictx->ctx->gdist->F(kA2, ictx->Ya);
    double dotfac = (ictx->q1x * ictx->q2x + ictx->q1y * ictx->q2y) / (ictx->q12 * ictx->q22);
    double value = ictx->ctx->Nc * M_1_PI * ictx->alphas_2pi * ictx->ctx->Sperp * ictx->qqfactor / ictx->z2 * (1 + ictx->xi2) / ictx->xi
        * dotfac * Fg1 * ictx->Fkq2;
    checkfinite(value);
    *real = value;
    *imag = 0;
}

void H14qqDelta::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double kA02 = gsl_pow_2(ictx->q1x - (1 - ictx->xiprime) * ictx->kT) + gsl_pow_2(ictx->q1y); // (q - (1 - xi') k)^2
    double logfac0 = log(kA02 / ictx->kT2);
    double logfac1 = log(ictx->q12 / ictx->kT2);
    double value = -ictx->ctx->Nc * ictx->alphas_2pi * ictx->ctx->Sperp * ictx->qqfactor / ictx->z2
        * ((1 + ictx->xiprime2) * logfac0 - 2 * logfac1) / (1 - ictx->xiprime)
        * ictx->Fk * ictx->Fkq1;
    checkfinite(value);
    *real = value;
    *imag = 0;
}



void H1qqCorrection::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double Iq = -61./12. + 3 * ictx->xp + 1.5 * gsl_pow_2(ictx->xp) + gsl_pow_3(ictx->xp) / 3. + 0.25 * gsl_pow_4(ictx->xp) + 4 * log(1 - ictx->xp);
    double resummation_factor = ictx->ctx->resummation_constant * ictx->alphas_2pi * ictx->ctx->Nc * 0.5 * Iq;
    double value = ictx->qqfactor / ictx->z2 * ictx->q12 * ictx->Fq1 / (M_PI * gsl_pow_2(ictx->kT2))
        * (exp(resummation_factor) - 1 - resummation_factor);
    checkfinite(value);
    *real = value;
    *imag = 0;
}

void H1qqCorrectionA::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double Iq = -61./12. + 3 * ictx->xp + 1.5 * gsl_pow_2(ictx->xp) + gsl_pow_3(ictx->xp) / 3. + 0.25 * gsl_pow_4(ictx->xp) + 4 * log(1 - ictx->xp);
    double resummation_factor = ictx->ctx->resummation_constant * ictx->alphas_2pi * ictx->ctx->Nc * 0.5 * Iq;
    double value = ictx->qqfactor / ictx->z2 * ictx->q12 * ictx->Fq1 / (M_PI * gsl_pow_2(ictx->kT2))
        * (-resummation_factor);
    checkfinite(value);
    *real = value;
    *imag = 0;
}

void H1qqCorrectionB::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double Iq = -61./12. + 3 * ictx->xp + 1.5 * gsl_pow_2(ictx->xp) + gsl_pow_3(ictx->xp) / 3. + 0.25 * gsl_pow_4(ictx->xp) + 4 * log(1 - ictx->xp);
    double resummation_factor = ictx->ctx->resummation_constant * ictx->alphas_2pi * ictx->ctx->Nc * 0.5 * Iq;
    double value = ictx->qqfactor / ictx->z2 * ictx->q12 * ictx->Fq1 / (M_PI * gsl_pow_2(ictx->kT2))
        * (exp(resummation_factor) - 1);
    checkfinite(value);
    *real = value;
    *imag = 0;
}

void H1qqExact::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    double value = ictx->alphas_2pi * M_1_PI * ictx->qqfactor / ictx->z2 * (1 + ictx->xi2) / (1 - ictx->xi)
        * (ictx->ctx->CF * gsl_pow_2(1 - ictx->xi) + ictx->ctx->Nc * ictx->xi) / gsl_pow_2(ictx->kT2) * ictx->q12 * ictx->Fq1;
    checkfinite(value);
    *real = value;
    *imag = 0;
}


void H02gg::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double kA2 = gsl_pow_2(ictx->q1x - ictx->kT) + gsl_pow_2(ictx->q1y); // (q - k)^2
    double value = ictx->ctx->Sperp * ictx->ggfactor / ictx->z2 * ictx->Fq1 * ictx->Fkq1;
    checkfinite(value);
    *real = value;
    *imag = 0;
}


void H12qqbar::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double Pfac = 1 - 2 * ictx->xiprime + 2 * ictx->xiprime2;
    double kA2 = gsl_pow_2(ictx->q1x - ictx->xiprime * ictx->kT) + gsl_pow_2(ictx->q1y); // (q - xi' k)^2
    double logfac = log(kA2 / ictx->kT2);
    double value = -2 * ictx->ctx->Nf * ictx->ctx->TR * ictx->ctx->Sperp * ictx->alphas_2pi * ictx->ggfactor / ictx->z2
        * Pfac * ictx->Fq1 * ictx->Fkq1 * logfac;
    checkfinite(value);
    *real = value;
    *imag = 0;
}

void H16ggSingular::Fs(const IntegrationContext* ictx, double* real, double* imag) const {
    double Pfac = gsl_pow_2(1 - ictx->xi + ictx->xi2) / ictx->xi2;
    double dotfac = (ictx->q1x * ictx->q2x + ictx->q1y * ictx->q2y) / (ictx->q12 * ictx->q22);
    double kA2 = gsl_pow_2(ictx->kT + ictx->q1x + ictx->q3x) + gsl_pow_2(ictx->q1y + ictx->q3y); // (k + q1 + q3)^2
    double kB2 = gsl_pow_2(ictx->kT / ictx->xi + ictx->q2x + ictx->q3x) + gsl_pow_2(ictx->q2y + ictx->q3y); // (k / ξ + q2 + q3)^2
    double Fg1 = ictx->ctx->gdist->F(kA2, ictx->Ya);
    double Fg2 = ictx->ctx->gdist->F(kB2, ictx->Ya);
    double value = -4 * M_1_PI * ictx->alphas_2pi * ictx->ctx->Nc * ictx->ctx->Sperp * ictx->ggfactor / ictx->z2
         * Pfac * dotfac * Fg1 * Fg2 * ictx->Fq3;
    checkfinite(value);
    *real = value;
    *imag = 0;
}

void H16ggDelta::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double kA02 = gsl_pow_2(ictx->q2x - ictx->xiprime * ictx->kT) + gsl_pow_2(ictx->q2y); // (q2 - xi' k)^2
    double kA12 = gsl_pow_2(ictx->q2x - ictx->kT) + gsl_pow_2(ictx->q2y); // (q2 - k)^2
    double kB2 = gsl_pow_2(ictx->q1x + ictx->kT) + gsl_pow_2(ictx->q1y); // (q1 + k)^2
    double kC2 = gsl_pow_2(ictx->q1x + ictx->q2x) + gsl_pow_2(ictx->q1y + ictx->q2y); // (q1 + q2)^2
    double Fg2 = ictx->ctx->gdist->F(kB2, ictx->Ya);
    double Fg3 = ictx->ctx->gdist->F(kC2, ictx->Ya);
    double log0 = log(kA02 / ictx->kT2);
    double log1 = log(kA12 / ictx->kT2);
    double value = -4 * ictx->alphas_2pi * ictx->ctx->Nc * ictx->ctx->Sperp * ictx->ggfactor / ictx->z2 * ictx->Fq1 * Fg2 * Fg3
                      * ((ictx->xiprime * log0 - log1) / (1 - ictx->xiprime)
                         + 0.5 * ictx->xiprime * (1 - ictx->xiprime) * log0);
    checkfinite(value);
    *real = value;
    *imag = 0;
}

void H1ggCorrection::Fd(const IntegrationContext* ictx, double* real, double* imag) const {
    double Ig = -152./15. + 12 * ictx->xp - 6 * gsl_pow_2(ictx->xp) + 16 * gsl_pow_3(ictx->xp) / 3. - 2 * gsl_pow_4(ictx->xp) + 4 * gsl_pow_5(ictx->xp) / 5. + 8 * gsl_atanh(1 - 2 * ictx->xp);
    double resummation_factor = ictx->ctx->resummation_constant * ictx->alphas_2pi * ictx->ctx->Nc * Ig;
    double value = ictx->ggfactor / ictx->z2 * ictx->q12 * ictx->Fq1 / (M_PI * gsl_pow_2(ictx->kT2))
        * (exp(resummation_factor) - 1 - resummation_factor);
    checkfinite(value);
    *real = value;
    *imag = 0;
}

void H1ggExact::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    double value = ictx->alphas_2pi * ictx->ctx->Nc * M_1_PI * ictx->ggfactor / ictx->z2
        * 2 * gsl_pow_2(1 - ictx->xi + ictx->xi2) * (1 + ictx->xi2 + gsl_pow_2(1 - ictx->xi)) / (ictx->xi - ictx->xi2)
        / gsl_pow_2(ictx->kT2) * ictx->q12 * ictx->Fq1;
    checkfinite(value);
    *real = value;
    *imag = 0;
}

void H14gq::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    double Pfac = (2 - 2 * ictx->xi + ictx->xi2) / ictx->xi2; // Pgq(xi) / xi
    double dotfac = (ictx->q1x * ictx->q2x + ictx->q1y * ictx->q2y) / (ictx->q12 * ictx->q22);
    double kA2 = gsl_pow_2(ictx->kT / ictx->xi + ictx->q1x) + gsl_pow_2(ictx->q1y);
    double kB2 = gsl_pow_2(ictx->kT * (ictx->xi - 1) / ictx->xi - ictx->q1x + ictx->q2x) + gsl_pow_2(-ictx->q1y + ictx->q2y);
    double Fg1 = ictx->ctx->gdist->F(kA2, ictx->Ya);
    double Fg2 = ictx->ctx->gdist->F(kB2, ictx->Ya);
    double value = -ictx->alphas_2pi * ictx->ctx->Nc * M_1_PI * ictx->gqfactor / ictx->z2 * Pfac * Fg1 * Fg2 * dotfac;
    checkfinite(value);
    *real = value;
    *imag = 0;
}

void H14qg::Fn(const IntegrationContext* ictx, double* real, double* imag) const {
    double Pfac = (1 - 2 * ictx->xi + 2 * ictx->xi2) / ictx->xi; // Pqg(xi) / xi
    double dotfac = (ictx->q1x * ictx->q2x + ictx->q1y * ictx->q2y) / (ictx->q12 * ictx->q22);
    double kB2 = gsl_pow_2(ictx->kT * (ictx->xi - 1) / ictx->xi - ictx->q1x + ictx->q2x) + gsl_pow_2(-ictx->q1y + ictx->q2y);
    double Fg2 = ictx->ctx->gdist->F(kB2, ictx->Ya);
    double value = -ictx->alphas_2pi * M_1_PI * ictx->qgfactor / ictx->z2 * Pfac * ictx->Fkq1 * Fg2 * dotfac;
    checkfinite(value);
    *real = value;
    *imag = 0;
}
