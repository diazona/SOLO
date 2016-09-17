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
#include "integrationregion.h"

#define checkfinite(d) assert(gsl_finite(d))

static inline double zmax(const Context& ctx) {
    return 1;
}
/* EXACT_LIMIT_SCHEME == 0 means to use Bowen's set of integral limits for
 * exact kinematics; any other value or undefined means to use mine
 *
 * The explicit check for defined(EXACT_LIMIT_SCHEME) may be redundant but
 * I leave it here as a reminder in case the default changes
 */
static inline double zmin(const Context& ctx) {
#if defined(EXACT_LIMIT_SCHEME) && EXACT_LIMIT_SCHEME == 0
    return ctx.exact_kinematics ? ctx.tauhat : ctx.tau;
#else
    return ctx.tau;
#endif
}
static inline double ximax(const Context& ctx, const double z) {
#if defined(EXACT_LIMIT_SCHEME) && EXACT_LIMIT_SCHEME == 0
    const double xahat = sqrt(ctx.pT2) / (z * ctx.sqs) * exp(-ctx.Y);
    return ctx.exact_kinematics ? 1 - xahat : 1;
#else
    return 1;
#endif
}
static inline double ximin(const Context& ctx, const double z) {
    return ctx.tau / z;
}
static inline double ymax(const Context& ctx) {
    return 1;
}
static inline double ymin(const Context& ctx) {
    return ctx.tau;
}

static double linear_transform(double input, double input_min, double input_max, double output_min, double output_max) {
    // Special-case the min and max; sometimes they get messed up by floating-point roundoff error
    if (input == input_min) {
        return output_min;
    }
    else if (input == input_max) {
        return output_max;
    }
    else if (input_min == output_min && input_max == output_max) {
        return input;
    }
    else {
        return (input - input_min) * (output_max - output_min) / (input_max - input_min) + output_min;
    }
}

static double xi_zy(const Context& ctx, size_t core_dimensions, double z, double y) {
    if (core_dimensions == 1) {
        /* When calculating LO terms, or NLO terms in exact kinematics, this shouldn't matter
         * When calculating NLO terms using approximate kinematics, this needs to be 1
         */
        return 1;
    }
    else {
        assert(y <= 1);
        assert(y >= ctx.tau);
        double xi;
        if (y == ymax(ctx)) {
            // if y == 1 then the formula should set xi = 1 but sometimes it doesn't
            // because of floating point roundoff error, so do that manually
            xi = ximax(ctx, z);
        }
        else {
            xi = (y - ymin(ctx)) * (ximax(ctx, z) - ximin(ctx, z)) / (ymax(ctx) - ymin(ctx)) + ximin(ctx, z);
        }
        assert(ctx.exact_kinematics ? (xi < 1) : (xi <= 1));
        return xi;
    }
}

IntegrationRegion::IntegrationRegion(const size_t dimensions) : dimensions(dimensions) {}

IntegrationRegion::~IntegrationRegion() {}


/*
 * In general, the values passed as the third argument to update() are interpreted as
 *  values[0]: z
 *  values[1]: y (if applicable)
 *  values[2]: rx or sx or q1x
 *  values[3]: ry or sy or q1y
 *  values[4]: tx or q2x
 * and so on.
 */

// IntegrationRegion operator*(const IntegrationRegion& a, const IntegrationRegion& b) {
//     return CompositeIntegrationRegion(a, b);
// }

// default implementation
bool IntegrationRegion::operator<(const IntegrationRegion& other) const {
    // By default, all instances of a given IntegrationRegion subclass compare equal
    return typeid(*this).before(typeid(other));
}


size_t total_dimensions(const size_t n_subregions, const IntegrationRegion* const* subregions) {
    size_t dimensions = 0;
    for (size_t i = 0; i < n_subregions; i++) {
        dimensions += subregions[i]->dimensions;
    }
    return dimensions;
}

// This allows having more than two constructor arguments if desired
vector<const IntegrationRegion*> combine_integration_regions(size_t n, const IntegrationRegion& a, const IntegrationRegion& b) {
    // TODO: if an argument is a CompositeIntegrationRegion, add its components instead of the region itself
    vector<const IntegrationRegion*> vec;
    switch (n) { // intentionally omitting break statements
        case 1:
            vec.push_back(&a);
        case 2:
            vec.push_back(&b);
        // easy to continue the pattern to more arguments
        default:
            assert(false);
    }
    return vec;
}

CompositeIntegrationRegion::CompositeIntegrationRegion() :
  IntegrationRegion(0) {
}

CompositeIntegrationRegion::CompositeIntegrationRegion(const IntegrationRegion& a, const IntegrationRegion& b) :
  IntegrationRegion(a.dimensions + b.dimensions),
  m_subregions(combine_integration_regions(2, a, b)) {
}

CompositeIntegrationRegion::CompositeIntegrationRegion(const vector<const IntegrationRegion*> subregions) :
  IntegrationRegion(total_dimensions(subregions.size(), &subregions[0])),
  m_subregions(subregions) {
}

CompositeIntegrationRegion::CompositeIntegrationRegion(const size_t n_subregions, const IntegrationRegion** subregions) :
  IntegrationRegion(total_dimensions(n_subregions, subregions)),
  m_subregions(subregions, subregions + n_subregions) {
}

void CompositeIntegrationRegion::fill_min(const Context& ctx, double* min) const {
    for (size_t k = 0; k < m_subregions.size(); k++) {
        m_subregions[k]->fill_min(ctx, min);
        min += m_subregions[k]->dimensions;
    }
}

void CompositeIntegrationRegion::fill_max(const Context& ctx, double* max) const {
    for (size_t k = 0; k < m_subregions.size(); k++) {
        m_subregions[k]->fill_max(ctx, max);
        max += m_subregions[k]->dimensions;
    }
}

double CompositeIntegrationRegion::jacobian(const IntegrationContext &ictx) const {
    double jacobian = 1;
    for (size_t k = 0; k < m_subregions.size(); k++) {
        jacobian *= m_subregions[k]->jacobian(ictx);
    }
    return jacobian;
}

void CompositeIntegrationRegion::update(IntegrationContext &ictx, const double *values) const {
    for (size_t k = 0; k < m_subregions.size(); k++) {
        m_subregions[k]->update(ictx, values);
        values += m_subregions[k]->dimensions;
    }
}

bool CompositeIntegrationRegion::operator<(const IntegrationRegion& other) const {
    const CompositeIntegrationRegion* p_other = dynamic_cast<const CompositeIntegrationRegion*>(&other);
    if (p_other != NULL) {
        if (m_subregions.size() < p_other->m_subregions.size()) {
            return true;
        }
        else if (m_subregions.size() == p_other->m_subregions.size()) {
            for (size_t i = 0; i < m_subregions.size(); i++) {
                const IntegrationRegion* i1 = m_subregions[i];
                const IntegrationRegion* i2 = p_other->m_subregions[i];
                if (*i1 < *i2) {
                    return true;
                }
                else if (*i2 < *i1) {
                    return false;
                }
            }
        }
        return false;
    }
    else {
        return IntegrationRegion::operator<(other);
    }
}

CompositeIntegrationRegion CompositeIntegrationRegion::operator*=(const IntegrationRegion& other) {
    return CompositeIntegrationRegion(*this, other);
}


LOKinematicsIntegrationRegion::LOKinematicsIntegrationRegion() :
  IntegrationRegion(1) {
}

void LOKinematicsIntegrationRegion::fill_min(const Context& ctx, double* min) const {
    min[0] = ctx.tau;
}

void LOKinematicsIntegrationRegion::fill_max(const Context& ctx, double* max) const {
    max[0] = 1;
}

double LOKinematicsIntegrationRegion::jacobian(const IntegrationContext& ictx) const {
    return 1;
}

void LOKinematicsIntegrationRegion::update(IntegrationContext& ictx, const double* values) const {
    ictx.z = values[0];
    ictx.xi = 1;
}


NLOKinematicsIntegrationRegion::NLOKinematicsIntegrationRegion() :
  IntegrationRegion(2) {
}

void NLOKinematicsIntegrationRegion::fill_min(const Context& ctx, double* min) const {
    min[0] = ctx.tau;
    min[1] = ctx.tau;
}

void NLOKinematicsIntegrationRegion::fill_max(const Context& ctx, double* max) const {
    max[0] = 1;
    max[1] = 1;
}

double NLOKinematicsIntegrationRegion::jacobian(const IntegrationContext& ictx) const {
    double jacobian = (1 - ictx.ctx.tau / ictx.z) / (1 - ictx.ctx.tau);
    checkfinite(jacobian);
    return jacobian;
}

void NLOKinematicsIntegrationRegion::update(IntegrationContext& ictx, const double* values) const {
    ictx.z = values[0];
    ictx.xi = linear_transform(values[1], ictx.ctx.tau, 1, ictx.ctx.tau / ictx.z, 1);
}


NLOClippedKinematicsIntegrationRegion::NLOClippedKinematicsIntegrationRegion() :
  IntegrationRegion(2) {
}

void NLOClippedKinematicsIntegrationRegion::fill_min(const Context& ctx, double* min) const {
    min[0] = ctx.tauhat;
    min[1] = ctx.tau;
}

void NLOClippedKinematicsIntegrationRegion::fill_max(const Context& ctx, double* max) const {
    max[0] = 1;
    max[1] = 1;
}

static inline double xahat(const IntegrationContext& ictx) {
    // This happens to be the same as xg, but I write it out explicitly here
    // just in case the definition of xg ever changes or something
    return sqrt(ictx.ctx.pT2) / (ictx.z * ictx.ctx.sqs) * exp(-ictx.ctx.Y);
}

double NLOClippedKinematicsIntegrationRegion::jacobian(const IntegrationContext& ictx) const {
    double jacobian = (1 - xahat(ictx) - ictx.ctx.tau / ictx.z) / (1 - ictx.ctx.tau);
    checkfinite(jacobian);
    return jacobian;
}

void NLOClippedKinematicsIntegrationRegion::update(IntegrationContext& ictx, const double* values) const {
    ictx.z = values[0];
    ictx.xi = linear_transform(values[1], ictx.ctx.tau, 1, ictx.ctx.tau / ictx.z, 1 - xahat(ictx));
}


CartesianIntegrationRegion::CartesianIntegrationRegion() :
  IntegrationRegion(2) {
}

void CartesianIntegrationRegion::fill_min(const Context& ctx, double* min) const {
    min[0] = -ctx.inf;
    min[1] = -ctx.inf;
}

void CartesianIntegrationRegion::fill_max(const Context& ctx, double* max) const {
    max[0] = ctx.inf;
    max[1] = ctx.inf;
}

double CartesianIntegrationRegion::jacobian(const IntegrationContext& ictx) const {
    return 1;
}

void CartesianIntegrationRegion::update(IntegrationContext& ictx, const double* values) const {
    x(ictx) = values[0];
    y(ictx) = values[1];
}


PolarIntegrationRegion::PolarIntegrationRegion() :
  IntegrationRegion(2) {
}


PolarIntegrationRegion::PolarIntegrationRegion(const size_t dimensions) :
  IntegrationRegion(dimensions) {
}

void PolarIntegrationRegion::fill_min(const Context& ctx, double* min) const {
    min[0] = ctx.cutoff;
    min[1] = 0;
}

void PolarIntegrationRegion::fill_max(const Context& ctx, double* max) const {
    max[0] = ctx.inf;
    max[1] = 2 * M_PI;
}

double PolarIntegrationRegion::jacobian(const IntegrationContext& ictx) const {
    double jacobian = r(ictx);
    checkfinite(jacobian);
    return jacobian;
}

void PolarIntegrationRegion::update(IntegrationContext& ictx, const double* values) const {
    double r = values[0];
    double theta = values[1];
    double center_point[2];
    center(ictx, center_point);
    x(ictx) = r * cos(theta) + center_point[0];
    y(ictx) = r * sin(theta) + center_point[1];
}

// default implementation
double PolarIntegrationRegion::r(const IntegrationContext& ictx) const {
    return gsl_hypot(x(ictx), y(ictx));
}

// default implementation
void PolarIntegrationRegion::center(const IntegrationContext& ictx, double* center_point) const {
    center_point[0] = center_point[1] = 0;
}


RadialIntegrationRegion::RadialIntegrationRegion() : PolarIntegrationRegion() {}

void RadialIntegrationRegion::fill_min(const Context& ctx, double* min) const {
    min[0] = 0;
}

void RadialIntegrationRegion::fill_max(const Context& ctx, double* max) const {
    max[0] = ctx.inf;
}

double RadialIntegrationRegion::jacobian(const IntegrationContext& ictx) const {
    return PolarIntegrationRegion::jacobian(ictx) * 2 * M_PI;
}

void RadialIntegrationRegion::update(IntegrationContext& ictx, const double* values) const {
    double r_and_theta[2];
    r_and_theta[0] = values[0];
    r_and_theta[1] = 0;
    PolarIntegrationRegion::update(ictx, r_and_theta);
}


void PolarLimitedIntegrationRegion::fill_max(const Context& ctx, double* max) const {
    max[0] = 1;
    max[1] = 2 * M_PI;
}

double PolarLimitedIntegrationRegion::jacobian(const IntegrationContext& ictx) const {
    double jacobian = PolarIntegrationRegion::jacobian(ictx) * limit(ictx);
    checkfinite(jacobian);
    return jacobian;
}

void PolarLimitedIntegrationRegion::update(IntegrationContext& ictx, const double* values) const {
    // these should be consistent with the formulas in IntegrationContext::recalculate_everything()
    ictx.kT2 = ictx.ctx.pT2 / ictx.z2;
    ictx.kT = sqrt(ictx.kT2);
    ictx.xp = ictx.ctx.tau / (ictx.z * ictx.xi);

    double r_and_theta[2];
    r_and_theta[0] = values[0] * limit(ictx);
    r_and_theta[1] = values[1];
    PolarIntegrationRegion::update(ictx, r_and_theta);
}


double ExactKinematicIntegrationRegion::limit(const IntegrationContext& ictx) const {
    double lim1 = sqrt(ictx.kT * (ictx.ctx.sqs * exp(ictx.ctx.Y) - ictx.kT) * (1 - ictx.xi) / ictx.xi);
    // quick and dirty way to check the equivalence of these two formulas while testing
#ifndef NDEBUG
    if (xi != 1) {
        double lim2 = sqrt((xp * ictx.ctx.sqs * ictx.ctx.sqs - kT2) * (1 - xi) / xi);
        assert(gsl_fcmp(lim1, lim2, 1e-10) == 0);
    }
#endif
    return lim1;
}

void ExactKinematicIntegrationRegion::center(const IntegrationContext& ictx, double* center_point) const {
    center_point[0] = ictx.kT;
    center_point[1] = 0;
}


/*
double& R1IntegrationRegion::x(IntegrationContext& ictx) const {
    return ictx.xx;
}

double& R1IntegrationRegion::y(IntegrationContext& ictx) const {
    return ictx.xy;
}

double& R2IntegrationRegion::x(IntegrationContext& ictx) const {
    return ictx.yx;
}

double& R2IntegrationRegion::y(IntegrationContext& ictx) const {
    return ictx.yy;
}

double& R3IntegrationRegion::x(IntegrationContext& ictx) const {
    return ictx.bx;
}

double& R3IntegrationRegion::y(IntegrationContext& ictx) const {
    return ictx.by;
}

double& Q1IntegrationRegion::x(IntegrationContext& ictx) const {
    return ictx.q1x;
}

double& Q1IntegrationRegion::y(IntegrationContext& ictx) const {
    return ictx.q1y;
}

double& Q2IntegrationRegion::x(IntegrationContext& ictx) const {
    return ictx.q2x;
}

double& Q2IntegrationRegion::y(IntegrationContext& ictx) const {
    return ictx.q2y;
}

double& Q3IntegrationRegion::x(IntegrationContext& ictx) const {
    return ictx.q3x;
}

double& Q3IntegrationRegion::y(IntegrationContext& ictx) const {
    return ictx.q3y;
}

double& R1PIntegrationRegion::x(IntegrationContext& ictx) const {
    return ictx.xx;
}

double& R1PIntegrationRegion::y(IntegrationContext& ictx) const {
    return ictx.xy;
}

double& R2PIntegrationRegion::x(IntegrationContext& ictx) const {
    return ictx.yx;
}

double& R2PIntegrationRegion::y(IntegrationContext& ictx) const {
    return ictx.yy;
}

double& R3PIntegrationRegion::x(IntegrationContext& ictx) const {
    return ictx.bx;
}

double& R3PIntegrationRegion::y(IntegrationContext& ictx) const {
    return ictx.by;
}

double& Q1PIntegrationRegion::x(IntegrationContext& ictx) const {
    return ictx.q1x;
}

double& Q1PIntegrationRegion::y(IntegrationContext& ictx) const {
    return ictx.q1y;
}

const double Q1PIntegrationRegion::r(const IntegrationContext& ictx) const {
    return sqrt(ictx.q12);
}

double& Q2PIntegrationRegion::x(IntegrationContext& ictx) const {
    return ictx.q2x;
}

double& Q2PIntegrationRegion::y(IntegrationContext& ictx) const {
    return ictx.q2y;
}

const double Q2PIntegrationRegion::r(const IntegrationContext& ictx) const {
    return sqrt(ictx.q22);
}

double& Q3PIntegrationRegion::x(IntegrationContext& ictx) const {
    return ictx.q3x;
}

double& Q3PIntegrationRegion::y(IntegrationContext& ictx) const {
    return ictx.q3y;
}

const double Q3PIntegrationRegion::r(const IntegrationContext& ictx) const {
    return sqrt(ictx.q32);
}*/

// void AngleIndependentPositionIntegrationRegion::fill_min(const Context& ctx, const size_t core_dimensions, double* min) const {
//     assert(ctx.tau < 1);
//     assert(ctx.tauhat < 1);
//     assert(core_dimensions == 1 || core_dimensions == 2);
//     size_t i = 0;
//     min[i++] = zmin(ctx);
//     if (core_dimensions == 2) {
//         min[i++] = ymin(ctx);
//     }
//     while (i < core_dimensions + dimensions) { min[i++] = 0; }
// }
//
// void AngleIndependentPositionIntegrationRegion::fill_max(const Context& ctx, const size_t core_dimensions, double* max) const {
//     assert(core_dimensions == 1 || core_dimensions == 2);
//     size_t i = 0;
//     max[i++] = zmax(ctx);
//     if (core_dimensions == 2) {
//         max[i++] = ymax(ctx);
//     }
//     while (i < core_dimensions + dimensions) {
//         max[i++] = ctx.inf;
//     }
// }
//
// double AngleIndependentPositionIntegrationRegion::jacobian(const unsigned int ncoords, const double* coordinates, const IntegrationContext& ictx, const size_t core_dimensions) const {
//     assert(dimensions <= 3);
//     double jacobian = IntegrationRegion::jacobian(ncoords, coordinates, ictx, core_dimensions); // y to xi
//     // now (r,theta) to (x,y)
//     if (dimensions > 0) {
//         assert(ictx.xy == 0);
//         jacobian *= 2 * M_PI * ictx.xx;
//     }
//     if (dimensions > 1) {
//         assert(ictx.yy == 0);
//         jacobian *= 2 * M_PI * ictx.yx;
//     }
//     if (dimensions > 2) {
//         assert(ictx.by == 0);
//         jacobian *= 2 * M_PI * ictx.bx;
//     }
//     checkfinite(jacobian);
//     return jacobian;
// }
// void AngleIndependentPositionIntegrationRegion::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
//     assert(core_dimensions == 1 || core_dimensions == 2);
//     ictx.update_kinematics(values[0], xi_zy(ictx.ctx, core_dimensions, values[0], values[1]), core_dimensions);
//     ictx.update_positions(
//         dimensions > 0 ? values[core_dimensions + 0] : 0,
//         0,
//         dimensions > 1 ? values[core_dimensions + 1] : 0,
//         0,
//         dimensions > 2 ? values[core_dimensions + 2] : 0,
//         0
//     );
//     ictx.update_parton_functions();
// }
//
// double RescaledAngleIndependentPositionIntegrationRegion::jacobian(const unsigned int ncoords, const double* coordinates, const IntegrationContext& ictx, const size_t core_dimensions) const {
//     assert(dimensions <= 3);
//     double jacobian = IntegrationRegion::jacobian(ncoords, coordinates, ictx, core_dimensions); // y to xi
//     // now (r,theta) to (x,y)
//     if (dimensions > 0) {
//         assert(ictx.xy == 0);
//         jacobian *= 2 * M_PI * ictx.xx * ictx.xi;
//     }
//     if (dimensions > 1) {
//         assert(ictx.yy == 0);
//         jacobian *= 2 * M_PI * ictx.yx * ictx.xi;
//     }
//     if (dimensions > 2) {
//         assert(ictx.by == 0);
//         jacobian *= 2 * M_PI * ictx.bx * ictx.xi;
//     }
//     checkfinite(jacobian);
//     return jacobian;
// }
// void RescaledAngleIndependentPositionIntegrationRegion::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
//     assert(core_dimensions == 1 || core_dimensions == 2);
//     ictx.update_kinematics(values[0], xi_zy(ictx.ctx, core_dimensions, values[0], values[1]), core_dimensions);
//     ictx.update_positions(
//         dimensions > 0 ? values[core_dimensions + 0] * ictx.xi : 0,
//         0,
//         dimensions > 1 ? values[core_dimensions + 1] * ictx.xi : 0,
//         0,
//         dimensions > 2 ? values[core_dimensions + 2] * ictx.xi : 0,
//         0
//     );
//     ictx.update_parton_functions();
// }
//
// void QLimitedMomentumIntegrationRegion::fill_max(const Context& ctx, const size_t core_dimensions, double* max) const {
//     assert(core_dimensions == 1 || core_dimensions == 2);
//     size_t i = 0;
//     max[i++] = zmax(ctx);
//     if (core_dimensions == 2) {
//         max[i++] = ymax(ctx);
//     }
//     while (i < core_dimensions + dimensions) {
//         max[i++] = 1;
//         max[i++] = 2 * M_PI;
//     }
// }
//
// double QLimitedMomentumIntegrationRegion::jacobian(const unsigned int ncoords, const double* coordinates, const IntegrationContext& ictx, const size_t core_dimensions) const {
//     double jacobian = RadialMomentumIntegrationRegion::jacobian(ncoords, coordinates, ictx, core_dimensions);
//     for (size_t i = 0; i < dimensions; i += 2) {
//         jacobian *= ictx.qmax;
//     }
//     checkfinite(jacobian);
//     return jacobian;
// }
//
// void QLimitedMomentumIntegrationRegion::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
//     // the thing in the array is a scaled integration variable between 0 and 1, so convert it to q
//     assert(core_dimensions == 1 || core_dimensions == 2);
//     ictx.update_kinematics(values[0], xi_zy(ictx.ctx, core_dimensions, values[0], values[1]), core_dimensions);
//     // qmax is computed in update_kinematics
//     double qmax = ictx.qmax;
//     ictx.update_momenta(
//         dimensions > 0 ? qmax * values[core_dimensions + 0] * cos(values[core_dimensions + 1]) + ictx.kT : 0,
//         dimensions > 0 ? qmax * values[core_dimensions + 0] * sin(values[core_dimensions + 1]) : 0,
//         dimensions > 2 ? qmax * values[core_dimensions + 2] * cos(values[core_dimensions + 3]) + ictx.kT : 0,
//         dimensions > 2 ? qmax * values[core_dimensions + 2] * sin(values[core_dimensions + 3]) : 0,
//         dimensions > 4 ? qmax * values[core_dimensions + 4] * cos(values[core_dimensions + 5]) + ictx.kT : 0,
//         dimensions > 4 ? qmax * values[core_dimensions + 4] * sin(values[core_dimensions + 5]) : 0
//     );
//     ictx.update_parton_functions();
// }
