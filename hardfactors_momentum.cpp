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
    *real = value;
    *imag = 0;
}

void H16ggSingular::Fs(IntegrationContext* ictx, double* real, double* imag) {
    double value = -4 * M_1_PI * ictx->alphasbar * ictx->ctx->Nc * ictx->ctx->Sperp * ictx->ggfactor / ictx->z2 *
        (1 - ictx->xi + ictx->xi2)*(1 - ictx->xi + ictx->xi2) / ictx->xi2 *
        (ictx->q1x * ictx->q2x + ictx->q1y * ictx->q2y) / (ictx->q12 * ictx->q22) *
        ictx->ctx->gdist->F((ictx->kT + ictx->q1x + ictx->q3x)*(ictx->kT + ictx->q1x + ictx->q3x) + (ictx->q1y + ictx->q3y)*(ictx->q1y + ictx->q3y), ictx->Qs2) *
        ictx->ctx->gdist->F((ictx->q1x / ictx->xi + ictx->q2x + ictx->q3x)*(ictx->q1x / ictx->xi + ictx->q2x + ictx->q3x) + (ictx->q1y / ictx->xi + ictx->q2y + ictx->q3y)*(ictx->q1y / ictx->xi + ictx->q2y + ictx->q3y), ictx->Qs2) *
        ictx->Fq3;
    *real = value;
    *imag = 0;
}

void H16ggDelta::Fd(IntegrationContext* ictx, double* real, double* imag) {
    double value = -4 * ictx->alphasbar * ictx->ctx->Nc * ictx->ctx->Sperp * ictx->ggfactor / ictx->z2 * (
        ( ictx->xi * log(((ictx->q2x - ictx->xi * ictx->kT)*(ictx->q2x - ictx->xi * ictx->kT) + ictx->q2y*ictx->q2y) / ictx->kT2) -
                     log(((ictx->q2x - ictx->kT)*(ictx->q2x - ictx->kT) + ictx->q2y*ictx->q2y) / ictx->kT2) ) / (1 - ictx->xi)
          + 0.5 * ictx->xi * (1 - ictx->xi) * log(((ictx->q2x - ictx->xi * ictx->kT)*(ictx->q2x - ictx->xi * ictx->kT) + ictx->q2y*ictx->q2y) / ictx->kT2)
        ) * ictx->Fq1
          * ictx->ctx->gdist->F((ictx->kT + ictx->q1x)*(ictx->kT + ictx->q1x) + ictx->q1y * ictx->q1y, ictx->Qs2)
          * ictx->ctx->gdist->F((ictx->q1x + ictx->q2x)*(ictx->q1x + ictx->q2x) + (ictx->q1y + ictx->q2y)*(ictx->q1y + ictx->q2y), ictx->Qs2);
    *real = value;
    *imag = 0;
}

void H14gq::Fn(IntegrationContext* ictx, double* real, double* imag) {
    double value = -ictx->alphasbar * ictx->ctx->Nc * M_1_PI * (2 - 2 * ictx->xi + ictx->xi2) / ictx->xi2 *
        ictx->ctx->gdist->F((ictx->kT / ictx->xi + ictx->q1x)*(ictx->kT / ictx->xi + ictx->q1x) + ictx->q1y*ictx->q1y, ictx->Qs2) *
        ictx->ctx->gdist->F((ictx->kT * (1 - 1. / ictx->xi) - ictx->q1x + ictx->q2x)*(ictx->kT * (1 - 1. / ictx->xi) - ictx->q1x + ictx->q2x) + (-ictx->q1y + ictx->q2y)*(-ictx->q1y + ictx->q2y), ictx->Qs2) *
        (ictx->q1x * ictx->q2x + ictx->q1y * ictx->q2y) / (ictx->q12 * ictx->q22);
    *real = value;
    *imag = 0;
}

void H14qg::Fn(IntegrationContext* ictx, double* real, double* imag) {
    double value = -ictx->alphasbar * M_1_PI * (1 - 2 * ictx->xi + 2 * ictx->xi2) / ictx->xi *
        ictx->ctx->gdist->F((ictx->kT - ictx->q1x)*(ictx->kT - ictx->q1x) + ictx->q1y*ictx->q1y, ictx->Qs2) *
        ictx->ctx->gdist->F((ictx->kT * (1 - 1. / ictx->xi) - ictx->q1x + ictx->q2x)*(ictx->kT * (1 - 1. / ictx->xi) - ictx->q1x + ictx->q2x) + (-ictx->q1y + ictx->q2y)*(-ictx->q1y + ictx->q2y), ictx->Qs2) *
        (ictx->q1x * ictx->q2x + ictx->q1y * ictx->q2y) / (ictx->q12 * ictx->q22);
    *real = value;
    *imag = 0;
}

