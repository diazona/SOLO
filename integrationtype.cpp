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

#include "integrationtype.h"

/**
 * For dimensions in which the lower and upper bounds are -infinity and +infinity,
 * we have to pick a finite value to cut off the integral. Using 10 seems to be
 * quite sufficient but this can always be raised. It shouldn't be made too large,
 * though, because if it is, the Monte Carlo sampling will miss the peak in the
 * integrand entirely and just output zero all the time.
 */
const double inf = 10;

/*
 * In general, the values passed as the third argument to update() are interpreted as
 *  values[0]: z
 *  values[1]: y (if applicable)
 *  values[2]: rx or sx or q1x
 *  values[3]: ry or sy or q1y
 *  values[4]: tx or q2x
 * and so on.
 */

void PlainIntegrationType::fill_min(IntegrationContext& ictx, const size_t core_dimensions, double* min) const {
    size_t i = 0;
    while (i < core_dimensions) { min[i++] = ictx.ctx->tau; }
    while (i < core_dimensions + extra_dimensions) { min[i++] = -inf; }
}
void PlainIntegrationType::fill_max(IntegrationContext& ictx, const size_t core_dimensions, double* max) const {
    size_t i = 0;
    while (i < core_dimensions) { max[i++] = 1; }
    while (i < core_dimensions + extra_dimensions) { max[i++] = inf; }
}

void DipoleIntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    if (core_dimensions == 1) {
        ictx.update_parton_factors(values[0], 1);
        ictx.update_positions(values[1], values[2], 0, 0, 0, 0);
    }
    else {
        ictx.update_parton_factors(values[0], values[1]);
        ictx.update_positions(values[2], values[3], 0, 0, 0, 0);
    }
}

void QuadrupoleIntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    if (core_dimensions == 1) {
        ictx.update_parton_factors(values[0], 1);
        ictx.update_positions(values[1], values[2], values[3], values[4], 0, 0);
    }
    else {
        ictx.update_parton_factors(values[0], values[1]);
        ictx.update_positions(values[2], values[3], values[4], values[5], 0, 0);
    }
}

void Momentum1IntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    if (core_dimensions == 1) {
        ictx.update_parton_factors(values[0], 1);
        ictx.update_momenta(values[1], values[2], 0, 0, 0, 0);
    }
    else {
        ictx.update_parton_factors(values[0], values[1]);
        ictx.update_momenta(values[2], values[3], 0, 0, 0, 0);
    }
}

void Momentum2IntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    if (core_dimensions == 1) {
        ictx.update_parton_factors(values[0], 1);
        ictx.update_momenta(values[1], values[2], values[3], values[4], 0, 0);
    }
    else {
        ictx.update_parton_factors(values[0], values[1]);
        ictx.update_momenta(values[2], values[3], values[4], values[5], 0, 0);
    }
}

void Momentum3IntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    if (core_dimensions == 1) {
        ictx.update_parton_factors(values[0], 1);
        ictx.update_momenta(values[1], values[2], values[3], values[4], values[5], values[6]);
    }
    else {
        ictx.update_parton_factors(values[0], values[1]);
        ictx.update_momenta(values[2], values[3], values[4], values[5], values[6], values[7]);
    }
}

/*
 * For these integration types, xiprime is taken to be the coordinate that
 * follows y (in the "2D" case) or z (in the "1D" case).
 */

void XiPIntegrationType::fill_min(IntegrationContext& ictx, const size_t core_dimensions, double* min) const {
    size_t i = 0;
    while (i < core_dimensions) { min[i++] = ictx.ctx->tau; }
    min[i++] = 0;
    while (i < core_dimensions + extra_dimensions) { min[i++] = -inf; }
}
void XiPIntegrationType::fill_max(IntegrationContext& ictx, const size_t core_dimensions, double* max) const {
    size_t i = 0;
    while (i < core_dimensions) { max[i++] = 1; }
    max[i++] = 1;
    while (i < core_dimensions + extra_dimensions) { max[i++] = inf; }
}

void Momentum1XiPIntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    if (core_dimensions == 1) {
        ictx.update_parton_factors(values[0], 1);
        ictx.update_auxiliary(values[1]);
        ictx.update_momenta(values[2], values[3], 0, 0, 0, 0);
    }
    else {
        ictx.update_parton_factors(values[0], values[1]);
        ictx.update_auxiliary(values[2]);
        ictx.update_momenta(values[3], values[4], 0, 0, 0, 0);
    }
}

void Momentum2XiPIntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    if (core_dimensions == 1) {
        ictx.update_parton_factors(values[0], 1);
        ictx.update_auxiliary(values[1]);
        ictx.update_momenta(values[2], values[3], values[4], values[5], 0, 0);
    }
    else {
        ictx.update_parton_factors(values[0], values[1]);
        ictx.update_auxiliary(values[2]);
        ictx.update_momenta(values[3], values[4], values[5], values[6], 0, 0);
    }
}
