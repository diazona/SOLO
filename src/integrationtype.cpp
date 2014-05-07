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
#include <typeinfo>
#include "integrationtype.h"

#define checkfinite(d) assert(gsl_finite(d))

bool compare_integration_types(const IntegrationType* a, const IntegrationType* b) {
    if (typeid(*a) == typeid(*b)) {
        return a->extra_dimensions < b->extra_dimensions;
    }
    else {
        return typeid(*a).before(typeid(*b));
    }
}

static inline double zmax(const Context* const ctx) {
    return 1;
}
/* EXACT_LIMIT_SCHEME == 0 means to use Bowen's set of integral limits for
 * exact kinematics; any other value or undefined means to use mine
 *
 * The explicit check for defined(EXACT_LIMIT_SCHEME) may be redundant but
 * I leave it here as a reminder in case the default changes
 */
static inline double zmin(const Context* const ctx) {
#if defined(EXACT_LIMIT_SCHEME) && EXACT_LIMIT_SCHEME == 0
    return ctx->exact_kinematics ? ctx->tauhat : ctx->tau;
#else
    return ctx->tau;
#endif
}
static inline double ximax(const Context* const ctx, const double z) {
#if defined(EXACT_LIMIT_SCHEME) && EXACT_LIMIT_SCHEME == 0
    const double xahat = sqrt(ctx->pT2) / (z * ctx->sqs) * exp(-ctx->Y);
    return ctx->exact_kinematics ? 1 - xahat : 1;
#else
    return 1;
#endif
}
static inline double ximin(const Context* const ctx, const double z) {
    return ctx->tau / z;
}
static inline double ymax(const Context* const ctx) {
    return 1;
}
static inline double ymin(const Context* const ctx) {
    return ctx->tau;
}

static double xi_zy(const Context* const ctx, size_t core_dimensions, double z, double y) {
    if (core_dimensions == 1) {
        /* When calculating LO terms, or NLO terms in exact kinematics, this shouldn't matter
         * When calculating NLO terms using approximate kinematics, this needs to be 1
         */
        return 1;
    }
    else {
        assert(y <= 1);
        assert(y >= ctx->tau);
        double xi;
        if (y == ymax(ctx)) {
            // if y == 1 then the formula should set xi = 1 but sometimes it doesn't
            // because of floating point roundoff error, so do that manually
            xi = ximax(ctx, z);
        }
        else {
            xi = (y - ymin(ctx)) * (ximax(ctx, z) - ximin(ctx, z)) / (ymax(ctx) - ymin(ctx)) + ximin(ctx, z);
        }
        assert(ctx->exact_kinematics ? (xi < 1) : (xi <= 1));
        return xi;
    }
}

/*
 * In general, the values passed as the third argument to update() are interpreted as
 *  values[0]: z
 *  values[1]: y (if applicable)
 *  values[2]: rx or sx or q1x
 *  values[3]: ry or sy or q1y
 *  values[4]: tx or q2x
 * and so on.
 */

// jacobian from y to xi
double IntegrationType::jacobian(IntegrationContext& ictx, const size_t core_dimensions) const {
    if (core_dimensions == 2) {
        double jacobian = (ximax(ictx.ctx, ictx.z) - ximin(ictx.ctx, ictx.z)) / (ymax(ictx.ctx) - ymin(ictx.ctx));
        checkfinite(jacobian);
        return jacobian;
    }
    else {
        return 1;
    }
}

void PlainIntegrationType::fill_min(const Context* const ctx, const size_t core_dimensions, double* min) const {
    assert(ctx->tau < 1);
    assert(ctx->tauhat < 1);
    assert(core_dimensions == 1 || core_dimensions == 2);
    size_t i = 0;
    min[i++] = zmin(ctx);
    if (core_dimensions == 2) {
        min[i++] = ymin(ctx);
    }
    while (i < core_dimensions + extra_dimensions) { min[i++] = -ctx->inf; }
}
void PlainIntegrationType::fill_max(const Context* const ctx, const size_t core_dimensions, double* max) const {
    assert(core_dimensions == 1 || core_dimensions == 2);
    size_t i = 0;
    max[i++] = zmax(ctx);
    if (core_dimensions == 2) {
        max[i++] = ymax(ctx);
    }
    while (i < core_dimensions + extra_dimensions) { max[i++] = ctx->inf; }
}

void PositionIntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    assert(core_dimensions == 1 || core_dimensions == 2);
    ictx.update_kinematics(values[0], xi_zy(ictx.ctx, core_dimensions, values[0], values[1]), core_dimensions);
    ictx.update_positions(
        extra_dimensions > 0 ? values[core_dimensions + 0] : 0,
        extra_dimensions > 1 ? values[core_dimensions + 1] : 0,
        extra_dimensions > 2 ? values[core_dimensions + 2] : 0,
        extra_dimensions > 3 ? values[core_dimensions + 3] : 0,
        extra_dimensions > 4 ? values[core_dimensions + 4] : 0,
        extra_dimensions > 5 ? values[core_dimensions + 5] : 0
    );
    ictx.update_parton_functions();
}
void MomentumIntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    assert(core_dimensions == 1 || core_dimensions == 2);
    ictx.update_kinematics(values[0], xi_zy(ictx.ctx, core_dimensions, values[0], values[1]), core_dimensions);
    ictx.update_momenta(
        extra_dimensions > 0 ? values[core_dimensions + 0] : 0,
        extra_dimensions > 1 ? values[core_dimensions + 1] : 0,
        extra_dimensions > 2 ? values[core_dimensions + 2] : 0,
        extra_dimensions > 3 ? values[core_dimensions + 3] : 0,
        extra_dimensions > 4 ? values[core_dimensions + 4] : 0,
        extra_dimensions > 5 ? values[core_dimensions + 5] : 0
    );
    ictx.update_parton_functions();
}

void NoIntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    ictx.update_kinematics(values[0], xi_zy(ictx.ctx, core_dimensions, values[0], values[1]), core_dimensions);
    ictx.update_parton_functions();
}



/*
 * For these integration types, xiprime is taken to be the coordinate that
 * follows y (in the "2D" case) or z (in the "1D" case).
 */

void XiPIntegrationType::fill_min(const Context* const ctx, const size_t core_dimensions, double* min) const {
    assert(ctx->tau < 1);
    assert(ctx->tauhat < 1);
    assert(core_dimensions == 1 || core_dimensions == 2);
    size_t i = 0;
    min[i++] = zmin(ctx);
    if (core_dimensions == 2) {
        min[i++] = ymin(ctx);
    }
    min[i++] = 0; // xiprime
    while (i < core_dimensions + extra_dimensions) { min[i++] = -ctx->inf; }
}
void XiPIntegrationType::fill_max(const Context* const ctx, const size_t core_dimensions, double* max) const {
    assert(core_dimensions == 1 || core_dimensions == 2);
    size_t i = 0;
    max[i++] = zmax(ctx);
    if (core_dimensions == 2) {
        max[i++] = ymax(ctx);
    }
    max[i++] = 1; // xiprime
    while (i < core_dimensions + extra_dimensions) { max[i++] = ctx->inf; }
}

void XiPIntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    assert(core_dimensions == 1 || core_dimensions == 2);
    ictx.update_kinematics(values[0], xi_zy(ictx.ctx, core_dimensions, values[0], values[1]), core_dimensions);
    ictx.update_auxiliary(values[core_dimensions]);
    ictx.update_momenta(
        extra_dimensions > 1 ? values[core_dimensions + 1] : 0,
        extra_dimensions > 2 ? values[core_dimensions + 2] : 0,
        extra_dimensions > 3 ? values[core_dimensions + 3] : 0,
        extra_dimensions > 4 ? values[core_dimensions + 4] : 0,
        extra_dimensions > 5 ? values[core_dimensions + 5] : 0,
        extra_dimensions > 6 ? values[core_dimensions + 6] : 0
    );
    ictx.update_parton_functions();
}

/*
 * For these integration types, values[2] and beyond are interpreted as magnitudes
 * of radius, with integration ranges 0 to infinity
 */
void RadialIntegrationType::fill_min(const Context* const ctx, const size_t core_dimensions, double* min) const {
    assert(ctx->tau < 1);
    assert(ctx->tauhat < 1);
    assert(core_dimensions == 1 || core_dimensions == 2);
    size_t i = 0;
    min[i++] = zmin(ctx);
    if (core_dimensions == 2) {
        min[i++] = ymin(ctx);
    }
    while (i < core_dimensions + extra_dimensions) { min[i++] = 0; }
}
void RadialIntegrationType::fill_max(const Context* const ctx, const size_t core_dimensions, double* max) const {
    assert(core_dimensions == 1 || core_dimensions == 2);
    size_t i = 0;
    max[i++] = zmax(ctx);
    if (core_dimensions == 2) {
        max[i++] = ymax(ctx);
    }
    while (i < core_dimensions + extra_dimensions) {
        max[i++] = ctx->inf;
        max[i++] = 2 * M_PI;
    }
}

double RadialPositionIntegrationType::jacobian(IntegrationContext& ictx, const size_t core_dimensions) const {
    assert(extra_dimensions <= 6);
    double jacobian = IntegrationType::jacobian(ictx, core_dimensions); // y to xi
    // now (r,theta) to (x,y)
    if (extra_dimensions > 0) {
        assert(extra_dimensions > 1);
        jacobian *= gsl_hypot(ictx.xx, ictx.xy);
    }
    if (extra_dimensions > 2) {
        assert(extra_dimensions > 3);
        jacobian *= gsl_hypot(ictx.yx, ictx.yy);
    }
    if (extra_dimensions > 4) {
        assert(extra_dimensions > 5);
        jacobian *= gsl_hypot(ictx.bx, ictx.by);
    }
    checkfinite(jacobian);
    return jacobian;
}
void RadialPositionIntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    assert(core_dimensions == 1 || core_dimensions == 2);
    ictx.update_kinematics(values[0], xi_zy(ictx.ctx, core_dimensions, values[0], values[1]), core_dimensions);
    ictx.update_positions(
        extra_dimensions > 0 ? values[core_dimensions + 0] * cos(values[core_dimensions + 1]) : 0,
        extra_dimensions > 0 ? values[core_dimensions + 0] * sin(values[core_dimensions + 1]) : 0,
        extra_dimensions > 2 ? values[core_dimensions + 2] * cos(values[core_dimensions + 3]) : 0,
        extra_dimensions > 2 ? values[core_dimensions + 2] * sin(values[core_dimensions + 3]) : 0,
        extra_dimensions > 4 ? values[core_dimensions + 4] * cos(values[core_dimensions + 5]) : 0,
        extra_dimensions > 4 ? values[core_dimensions + 4] * sin(values[core_dimensions + 5]) : 0
    );
    ictx.update_parton_functions();
}

double RadialMomentumIntegrationType::jacobian(IntegrationContext& ictx, const size_t core_dimensions) const {
    assert(extra_dimensions <= 6);
    double jacobian = IntegrationType::jacobian(ictx, core_dimensions); // y to xi
    // now (r,theta) to (x,y)
    if (extra_dimensions > 0) {
        assert(extra_dimensions > 1);
        jacobian *= sqrt(ictx.q12);
    }
    if (extra_dimensions > 2) {
        assert(extra_dimensions > 3);
        jacobian *= sqrt(ictx.q22);
    }
    if (extra_dimensions > 4) {
        assert(extra_dimensions > 5);
        jacobian *= sqrt(ictx.q32);
    }
    checkfinite(jacobian);
    return jacobian;
}
void RadialMomentumIntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    assert(core_dimensions == 1 || core_dimensions == 2);
    ictx.update_kinematics(values[0], xi_zy(ictx.ctx, core_dimensions, values[0], values[1]), core_dimensions);
    ictx.update_momenta(
        extra_dimensions > 0 ? values[core_dimensions + 0] * cos(values[core_dimensions + 1]) : 0,
        extra_dimensions > 0 ? values[core_dimensions + 0] * sin(values[core_dimensions + 1]) : 0,
        extra_dimensions > 2 ? values[core_dimensions + 2] * cos(values[core_dimensions + 3]) : 0,
        extra_dimensions > 2 ? values[core_dimensions + 2] * sin(values[core_dimensions + 3]) : 0,
        extra_dimensions > 4 ? values[core_dimensions + 4] * cos(values[core_dimensions + 5]) : 0,
        extra_dimensions > 4 ? values[core_dimensions + 4] * sin(values[core_dimensions + 5]) : 0
    );
    ictx.update_parton_functions();
}

void AngleIndependentPositionIntegrationType::fill_min(const Context*const ctx, const size_t core_dimensions, double* min) const {
    assert(ctx->tau < 1);
    assert(ctx->tauhat < 1);
    assert(core_dimensions == 1 || core_dimensions == 2);
    size_t i = 0;
    min[i++] = zmin(ctx);
    if (core_dimensions == 2) {
        min[i++] = ymin(ctx);
    }
    while (i < core_dimensions + extra_dimensions) { min[i++] = 0; }
}

void AngleIndependentPositionIntegrationType::fill_max(const Context* const ctx, const size_t core_dimensions, double* max) const {
    assert(core_dimensions == 1 || core_dimensions == 2);
    size_t i = 0;
    max[i++] = zmax(ctx);
    if (core_dimensions == 2) {
        max[i++] = ymax(ctx);
    }
    while (i < core_dimensions + extra_dimensions) {
        max[i++] = ctx->inf;
    }
}

double AngleIndependentPositionIntegrationType::jacobian(IntegrationContext& ictx, const size_t core_dimensions) const {
    assert(extra_dimensions <= 3);
    double jacobian = IntegrationType::jacobian(ictx, core_dimensions); // y to xi
    jacobian *= 2 * M_PI;
    // now (r,theta) to (x,y)
    if (extra_dimensions > 0) {
        assert(ictx.xy == 0);
        jacobian *= ictx.xx;
    }
    if (extra_dimensions > 1) {
        assert(ictx.yy == 0);
        jacobian *= ictx.yx;
    }
    if (extra_dimensions > 2) {
        assert(ictx.by == 0);
        jacobian *= ictx.bx;
    }
    checkfinite(jacobian);
    return jacobian;
}
void AngleIndependentPositionIntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    assert(core_dimensions == 1 || core_dimensions == 2);
    ictx.update_kinematics(values[0], xi_zy(ictx.ctx, core_dimensions, values[0], values[1]), core_dimensions);
    ictx.update_positions(
        extra_dimensions > 0 ? values[core_dimensions + 0] : 0,
        0,
        extra_dimensions > 1 ? values[core_dimensions + 1] : 0,
        0,
        extra_dimensions > 2 ? values[core_dimensions + 2] : 0,
        0
    );
    ictx.update_parton_functions();
}

void QLimitedMomentumIntegrationType::fill_max(const Context*const ctx, const size_t core_dimensions, double* max) const {
    assert(core_dimensions == 1 || core_dimensions == 2);
    size_t i = 0;
    max[i++] = zmax(ctx);
    if (core_dimensions == 2) {
        max[i++] = ymax(ctx);
    }
    while (i < core_dimensions + extra_dimensions) {
        max[i++] = 1;
        max[i++] = 2 * M_PI;
    }
}

double QLimitedMomentumIntegrationType::jacobian(IntegrationContext& ictx, const size_t core_dimensions) const {
    double jacobian = RadialMomentumIntegrationType::jacobian(ictx, core_dimensions);
    for (size_t i = 0; i < extra_dimensions; i += 2) {
        jacobian *= ictx.qmax;
    }
    checkfinite(jacobian);
    return jacobian;
}

void QLimitedMomentumIntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    // the thing in the array is a scaled integration variable between 0 and 1, so convert it to q
    assert(core_dimensions == 1 || core_dimensions == 2);
    ictx.update_kinematics(values[0], xi_zy(ictx.ctx, core_dimensions, values[0], values[1]), core_dimensions);
    // qmax is computed in update_kinematics
    double qmax = ictx.qmax;
    ictx.update_momenta(
        extra_dimensions > 0 ? qmax * values[core_dimensions + 0] * cos(values[core_dimensions + 1]) + ictx.kT : 0,
        extra_dimensions > 0 ? qmax * values[core_dimensions + 0] * sin(values[core_dimensions + 1]) : 0,
        extra_dimensions > 2 ? qmax * values[core_dimensions + 2] * cos(values[core_dimensions + 3]) + ictx.kT : 0,
        extra_dimensions > 2 ? qmax * values[core_dimensions + 2] * sin(values[core_dimensions + 3]) : 0,
        extra_dimensions > 4 ? qmax * values[core_dimensions + 4] * cos(values[core_dimensions + 5]) + ictx.kT : 0,
        extra_dimensions > 4 ? qmax * values[core_dimensions + 4] * sin(values[core_dimensions + 5]) : 0
    );
    ictx.update_parton_functions();
}
