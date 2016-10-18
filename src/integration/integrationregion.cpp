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


// IntegrationRegion::IntegrationRegion() {}
//
// IntegrationRegion::~IntegrationRegion() {}
//
// // default implementation
// bool IntegrationRegion::operator<(const IntegrationRegion& other) const {
//     // By default, all instances of a given IntegrationRegion subclass compare equal
//     return typeid(*this).before(typeid(other));
// }


// This allows having more constructor arguments if desired
template<class T>
vector<T*> make_vector(size_t n, T& a) {
    vector<T*> vec;
    switch (n) { // intentionally omitting break statements
        case 1:
            vec.push_back(&a);
        // easy to continue the pattern to more arguments
        default:
            assert(false);
    }
    return vec;
}

IntegrationRegion::IntegrationRegion(const CoreIntegrationRegion& core) :
  m_core_region(core) {
}

IntegrationRegion::IntegrationRegion(const CoreIntegrationRegion& core, const AuxiliaryIntegrationRegion& aux) :
  m_core_region(core),
  m_subregions(make_vector<const AuxiliaryIntegrationRegion>(1, aux)) {
}

IntegrationRegion::IntegrationRegion(const CoreIntegrationRegion& core, const vector<const AuxiliaryIntegrationRegion*> auxregions) :
  m_core_region(core),
  m_subregions(auxregions) {
}

IntegrationRegion::IntegrationRegion(const CoreIntegrationRegion& core, const size_t n_auxregions, const AuxiliaryIntegrationRegion** auxregions) :
  m_core_region(core),
  m_subregions(auxregions, auxregions + n_auxregions) {
}

IntegrationRegion::~IntegrationRegion() {}


void IntegrationRegion::fill_min(const Context& ctx, const bool xi_preintegrated_term, double* min) const {
    m_core_region.fill_min(ctx, xi_preintegrated_term, min);
    min += m_core_region.dimensions(xi_preintegrated_term);
    for (size_t k = 0; k < m_subregions.size(); k++) {
        m_subregions[k]->fill_min(ctx, min);
        min += m_subregions[k]->dimensions();
    }
}

void IntegrationRegion::fill_max(const Context& ctx, const bool xi_preintegrated_term, double* max) const {
    m_core_region.fill_max(ctx, xi_preintegrated_term, max);
    max += m_core_region.dimensions(xi_preintegrated_term);
    for (size_t k = 0; k < m_subregions.size(); k++) {
        m_subregions[k]->fill_max(ctx, max);
        max += m_subregions[k]->dimensions();
    }
}

double IntegrationRegion::jacobian(const IntegrationContext &ictx, const bool xi_preintegrated_term) const {
    double jacobian = m_core_region.jacobian(ictx, xi_preintegrated_term);
    for (size_t k = 0; k < m_subregions.size(); k++) {
        jacobian *= m_subregions[k]->jacobian(ictx);
    }
    return jacobian;
}

void IntegrationRegion::update(IntegrationContext &ictx, const bool xi_preintegrated_term, const double *values) const {
    m_core_region.update(ictx, xi_preintegrated_term, values);
    values += m_core_region.dimensions(xi_preintegrated_term);
    for (size_t k = 0; k < m_subregions.size(); k++) {
        m_subregions[k]->update(ictx, values);
        values += m_subregions[k]->dimensions();
    }
}

size_t IntegrationRegion::dimensions(const bool xi_preintegrated_term) const {
    size_t dim = m_core_region.dimensions(xi_preintegrated_term);
    for (size_t k = 0; k < m_subregions.size(); k++) {
        dim += m_subregions[k]->dimensions();
    }
    return dim;
}


bool IntegrationRegion::operator<(const IntegrationRegion& other) const {
    if (m_subregions.size() < other.m_subregions.size()) {
        return true;
    }
    else if (other.m_subregions.size() < m_subregions.size()) {
        return false;
    }
    else if (m_subregions.size() == other.m_subregions.size()) {
        if (m_core_region < other.m_core_region) {
            return true;
        }
        else if (other.m_core_region < m_core_region) {
            return false;
        }
        for (size_t i = 0; i < m_subregions.size(); i++) {
            const AuxiliaryIntegrationRegion* i1 = m_subregions[i];
            const AuxiliaryIntegrationRegion* i2 = other.m_subregions[i];
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

IntegrationRegion IntegrationRegion::operator*=(const AuxiliaryIntegrationRegion& other) {
    vector<const AuxiliaryIntegrationRegion*> auxregions(m_subregions);
    auxregions.push_back(&other);
    return IntegrationRegion(m_core_region, auxregions);
}


CoreIntegrationRegion::CoreIntegrationRegion() {}

CoreIntegrationRegion::~CoreIntegrationRegion() {}

bool CoreIntegrationRegion::operator<(const CoreIntegrationRegion& other) const {
    return typeid(*this).before(typeid(other));
}


AuxiliaryIntegrationRegion::AuxiliaryIntegrationRegion(const size_t dimensions) :
  m_dimensions(dimensions) {
}

AuxiliaryIntegrationRegion::~AuxiliaryIntegrationRegion() {}

size_t AuxiliaryIntegrationRegion::dimensions() const {
    return m_dimensions;
}

bool AuxiliaryIntegrationRegion::operator<(const AuxiliaryIntegrationRegion& other) const {
    return typeid(*this).before(typeid(other));
}


void LOKinematicsIntegrationRegion::fill_min(const Context& ctx, const bool xi_preintegrated_term, double* min) const {
    min[0] = ctx.tau;
}

void LOKinematicsIntegrationRegion::fill_max(const Context& ctx, const bool xi_preintegrated_term, double* max) const {
    max[0] = 1;
}

double LOKinematicsIntegrationRegion::jacobian(const IntegrationContext& ictx, const bool xi_preintegrated_term) const {
    return 1;
}

void LOKinematicsIntegrationRegion::update(IntegrationContext& ictx, const bool xi_preintegrated_term, const double* values) const {
    ictx.z = values[0];
    ictx.xi = 1;
}

size_t LOKinematicsIntegrationRegion::dimensions(const bool xi_preintegrated_term) const {
    return 1;
}

double LOKinematicsIntegrationRegion::effective_xi_min(const IntegrationContext& ictx) const {
    return 0;
}


NLOKinematicsIntegrationRegion::NLOKinematicsIntegrationRegion(const bool lower_zero) :
  m_lower_zero(lower_zero) {
}

void NLOKinematicsIntegrationRegion::fill_min(const Context& ctx, const bool xi_preintegrated_term, double* min) const {
    min[0] = ctx.tau;
    if (!xi_preintegrated_term) {
        min[1] = m_lower_zero ? 0 : ctx.tau;
    }
}

void NLOKinematicsIntegrationRegion::fill_max(const Context& ctx, const bool xi_preintegrated_term, double* max) const {
    max[0] = 1;
    if (!xi_preintegrated_term) {
        max[1] = 1;
    }
}

double NLOKinematicsIntegrationRegion::jacobian(const IntegrationContext& ictx, const bool xi_preintegrated_term) const {
    double jacobian = (xi_preintegrated_term || m_lower_zero) ? 1 : (1 - effective_xi_min(ictx)) / (1 - ictx.ctx.tau);
    checkfinite(jacobian);
    return jacobian;
}

void NLOKinematicsIntegrationRegion::update(IntegrationContext& ictx, const bool xi_preintegrated_term, const double* values) const {
    ictx.z = values[0];
    if (xi_preintegrated_term) {
        ictx.xi = 1;
    }
    else if (m_lower_zero) {
        ictx.xi = values[1];
    }
    else {
        ictx.xi = linear_transform(values[1], ictx.ctx.tau, 1, effective_xi_min(ictx), 1);
    }
}

size_t NLOKinematicsIntegrationRegion::dimensions(const bool xi_preintegrated_term) const {
    return xi_preintegrated_term ? 1 : 2;
}

double NLOKinematicsIntegrationRegion::effective_xi_min(const IntegrationContext& ictx) const {
    return m_lower_zero ? 0 : ictx.ctx.tau / ictx.z;
}


void NLOClippedKinematicsIntegrationRegion::fill_min(const Context& ctx, const bool xi_preintegrated_term, double* min) const {
    min[0] = ctx.tauhat;
    if (!xi_preintegrated_term) {
        min[1] = ctx.tau;
    }
}

void NLOClippedKinematicsIntegrationRegion::fill_max(const Context& ctx, const bool xi_preintegrated_term, double* max) const {
    max[0] = 1;
    if (!xi_preintegrated_term) {
        max[1] = 1;
    }
}

static inline double xahat(const IntegrationContext& ictx) {
    /* We have to use sqrt(pT2) / z instead of just using kT because this
     * function will be called from update(), which is called before
     * most of the values in the IntegrationContext are updated, and in
     * particular before kT is set to sqrt(pT2) / z. The only values in
     * ictx which we can count on to be updated are those which are set
     * right within update(), namely z.
     *
     * This formula also happens to be the same as that used for ictx.xg,
     * but that's basically a coincidence. It serves a different purpose
     * here.
     */
    return sqrt(ictx.ctx.pT2) / (ictx.z * ictx.ctx.sqs) * exp(-ictx.ctx.Y);
}

double NLOClippedKinematicsIntegrationRegion::jacobian(const IntegrationContext& ictx, const bool xi_preintegrated_term) const {
    double jacobian = xi_preintegrated_term ? 1 : (1 - xahat(ictx) - effective_xi_min(ictx)) / (1 - ictx.ctx.tau);
    checkfinite(jacobian);
    return jacobian;
}

void NLOClippedKinematicsIntegrationRegion::update(IntegrationContext& ictx, const bool xi_preintegrated_term, const double* values) const {
    ictx.z = values[0];
    ictx.xi = xi_preintegrated_term ? 1 : linear_transform(values[1], ictx.ctx.tau, 1, effective_xi_min(ictx), 1 - xahat(ictx));
}

size_t NLOClippedKinematicsIntegrationRegion::dimensions(const bool xi_preintegrated_term) const {
    return xi_preintegrated_term ? 1 : 2;
}

double NLOClippedKinematicsIntegrationRegion::effective_xi_min(const IntegrationContext& ictx) const {
    return ictx.ctx.tau / ictx.z;
}


CartesianIntegrationRegion::CartesianIntegrationRegion() :
  AuxiliaryIntegrationRegion(2) {
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
  AuxiliaryIntegrationRegion(2) {
}

PolarIntegrationRegion::PolarIntegrationRegion(const size_t dimensions) :
  AuxiliaryIntegrationRegion(dimensions) {
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


RadialIntegrationRegion::RadialIntegrationRegion() : PolarIntegrationRegion(1) {}

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
    double r_and_theta[2];
    r_and_theta[0] = values[0] * limit(ictx);
    r_and_theta[1] = values[1];
    PolarIntegrationRegion::update(ictx, r_and_theta);
}


double ExactKinematicIntegrationRegion::limit(const IntegrationContext& ictx) const {
    double z = ictx.z;
    double xi = ictx.xi;
    // these should be consistent with the formulas in IntegrationContext::recalculate_everything()
    double kT2 = ictx.ctx.pT2 / (z * z);
    double kT = sqrt(kT2);
    double xp = ictx.ctx.tau / z;
    double lim1 = sqrt(kT * (ictx.ctx.sqs * exp(ictx.ctx.Y) - kT) * (1 - xi) / xi);
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
