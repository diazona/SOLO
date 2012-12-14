// Hard factors in position space

#include <cassert>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include "hardfactor.h"
#include "hardfactors_momentum.h"

#define checkfinite(d) assert(gsl_finite(d))

using namespace momentum;

void H12qqbar::Fd(IntegrationContext* ictx, double* real, double* imag) {
    double value = 2 * ictx->ctx->Nf * ictx->ctx->TR * ictx->ctx->Sperp * ictx->alphasbar * ictx->ggfactor *
        (1 + 2 * ictx->xiprime + 2 * ictx->xiprime2) * ictx->Fq1 * ictx->Fkq1 * log(((ictx->q1x - ictx->xiprime * ictx->kT)*(ictx->q1x - ictx->xiprime * ictx->kT) + ictx->q1y*ictx->q1y)/ictx->kT2);
    checkfinite(value);
    *real = value;
    *imag = 0;
}

void H16ggSingular::Fs(IntegrationContext* ictx, double* real, double* imag) {
    double Pfac = gsl_pow_2(1 - ictx->xi + ictx->xi2) / ictx->xi2;
    double dotfac = (ictx->q1x * ictx->q2x + ictx->q1y * ictx->q2y) / (ictx->q12 * ictx->q22);
    double kA2 = gsl_pow_2(ictx->kT + ictx->q1x + ictx->q3x) + gsl_pow_2(ictx->q1y + ictx->q3y); // (k + q1 + q3)^2
    double kB2 = gsl_pow_2(ictx->kT / ictx->xi + ictx->q2x + ictx->q3x) + gsl_pow_2(ictx->q2y + ictx->q3y); // (k / ξ + q2 + q3)^2
    double Fg1 = ictx->ctx->gdist->F(kA2, ictx->Qs2);
    double Fg2 = ictx->ctx->gdist->F(kB2, ictx->Qs2);
    double Fg3 = ictx->ctx->gdist->F(ictx->q32, ictx->Qs2);
    double value = -4 * M_1_PI * ictx->alphasbar * ictx->ctx->Nc * ictx->ctx->Sperp * ictx->ggfactor / ictx->z2
         * Pfac * dotfac * Fg1 * Fg2 * Fg3;
    checkfinite(value);
    *real = value;
    *imag = 0;
}

void H16ggDelta::Fd(IntegrationContext* ictx, double* real, double* imag) {
    double kA02 = gsl_pow_2(ictx->q2x - ictx->xiprime * ictx->kT) + gsl_pow_2(ictx->q2y); // (q2 - xi' k)^2
    double kA12 = gsl_pow_2(ictx->q2x - ictx->kT) + gsl_pow_2(ictx->q2y); // (q2 - k)^2
    double kB2 = gsl_pow_2(ictx->q1x + ictx->kT) + gsl_pow_2(ictx->q1y); // (q1 + k)^2
    double kC2 = gsl_pow_2(ictx->q1x + ictx->q2x) + gsl_pow_2(ictx->q1y + ictx->q2y); // (q1 + q2)^2
    double Fg1 = ictx->ctx->gdist->F(ictx->q12, ictx->Qs2);
    double Fg2 = ictx->ctx->gdist->F(kB2, ictx->Qs2);
    double Fg3 = ictx->ctx->gdist->F(kC2, ictx->Qs2);
    double log0 = log(kA02 / ictx->kT2);
    double log1 = log(kA12 / ictx->kT2);
    double value = -4 * ictx->alphasbar * ictx->ctx->Nc * ictx->ctx->Sperp * ictx->ggfactor / ictx->z2 * Fg1 * Fg2 * Fg3
                      * ((ictx->xiprime * log0 - log1) / (1 - ictx->xiprime)
                         + 0.5 * ictx->xiprime * (1 - ictx->xiprime) * log0);
    checkfinite(value);
    *real = value;
    *imag = 0;
}

void H14gq::Fn(IntegrationContext* ictx, double* real, double* imag) {
    double value = -ictx->alphasbar * ictx->ctx->Nc * M_1_PI * (2 - 2 * ictx->xi + ictx->xi2) / ictx->xi2 *
        ictx->ctx->gdist->F((ictx->kT / ictx->xi + ictx->q1x)*(ictx->kT / ictx->xi + ictx->q1x) + ictx->q1y*ictx->q1y, ictx->Qs2) *
        ictx->ctx->gdist->F((ictx->kT * (1 - 1. / ictx->xi) - ictx->q1x + ictx->q2x)*(ictx->kT * (1 - 1. / ictx->xi) - ictx->q1x + ictx->q2x) + (-ictx->q1y + ictx->q2y)*(-ictx->q1y + ictx->q2y), ictx->Qs2) *
        (ictx->q1x * ictx->q2x + ictx->q1y * ictx->q2y) / (ictx->q12 * ictx->q22);
    checkfinite(value);
    *real = value;
    *imag = 0;
}

void H14qg::Fn(IntegrationContext* ictx, double* real, double* imag) {
    double value = -ictx->alphasbar * M_1_PI * (1 - 2 * ictx->xi + 2 * ictx->xi2) / ictx->xi *
        ictx->ctx->gdist->F((ictx->kT - ictx->q1x)*(ictx->kT - ictx->q1x) + ictx->q1y*ictx->q1y, ictx->Qs2) *
        ictx->ctx->gdist->F((ictx->kT * (1 - 1. / ictx->xi) - ictx->q1x + ictx->q2x)*(ictx->kT * (1 - 1. / ictx->xi) - ictx->q1x + ictx->q2x) + (-ictx->q1y + ictx->q2y)*(-ictx->q1y + ictx->q2y), ictx->Qs2) *
        (ictx->q1x * ictx->q2x + ictx->q1y * ictx->q2y) / (ictx->q12 * ictx->q22);
    checkfinite(value);
    *real = value;
    *imag = 0;
}

