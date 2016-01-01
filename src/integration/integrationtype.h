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

#ifndef _INTEGRATIONTYPE_H_
#define _INTEGRATIONTYPE_H_

#include "integrationcontext.h"

/**
 * Represents the properties of an integration region.
 *
 * An instance of this class is responsible for setting the number
 * of variables being integrated over and effectively specifying which
 * variables those are. In particular, it needs to set bounds on each
 * variable, and also translate the Monte Carlo integration variables
 * (which must have constant bounds) to the kinematic variables (whose
 * bounds can depend on each other) by calculating the jacobian determinant
 * and setting the value of the kinematic variables using provided values
 * of the integration variables.
 */
class IntegrationType {
public:
    /** The number of variables being integrated over excluding z and y */
    const size_t extra_dimensions;
    IntegrationType(const size_t extra_dimensions) : extra_dimensions(extra_dimensions) {}
    /**
     * Writes the lower bounds of the integration region into the array
     * `min` using information from the `Context`.
     */
    virtual void fill_min(const Context& ctx, const size_t core_dimensions, double* min) const = 0;
    /**
     * Writes the upper bounds of the integration region into the array
     * `max` using information from the `Context`.
     */
    virtual void fill_max(const Context& ctx, const size_t core_dimensions, double* max) const = 0;
    /**
     * Calculates the jacobian determinant to convert from the variables being
     * integrated over to the variables involved in the calculation.
     *
     * This gets called after `update()` at each step of the calculation.
     */
    virtual double jacobian(const unsigned int ncoords, const double* coordinates, const IntegrationContext& ictx, const size_t core_dimensions) const;
    /**
     * Updates the variables in the `IntegrationContext` based on the `values`,
     * which were passed in from the Monte Carlo integration routine.
     */
    virtual void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const = 0;
};

/**
 * Comparison function for IntegrationType pointers based on the
 * type_info and extra_dimensions fields.
 */
bool compare_integration_types(const IntegrationType* a, const IntegrationType* b);

/**
 * Base class for integration types which have one or two dimensions (z and/or y)
 * with limits of tau and 1, and any number of other dimensions with limits
 * plus or minus infinity
 */
class PlainIntegrationType : public IntegrationType {
protected:
    PlainIntegrationType(const size_t extra_dimensions) : IntegrationType(extra_dimensions) {}
    void fill_min(const Context& ctx, const size_t core_dimensions, double* min) const;
    void fill_max(const Context& ctx, const size_t core_dimensions, double* max) const;
};

/**
 * An integration type for integration over z and/or y only
 */
class NoIntegrationType : public PlainIntegrationType {
public:
    NoIntegrationType() : PlainIntegrationType(0) {}
    NoIntegrationType(NoIntegrationType const&);
private:
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

/**
 * Base class for Cartesian position integration types
 */
class PositionIntegrationType : public PlainIntegrationType {
public:
    PositionIntegrationType(const size_t extra_dimensions) : PlainIntegrationType(extra_dimensions) {}
protected:
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

/**
 * Base class for momentum integration types
 */
class MomentumIntegrationType : public PlainIntegrationType {
public:
    MomentumIntegrationType(const size_t extra_dimensions) : PlainIntegrationType(extra_dimensions) {}
protected:
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

/**
 * Base class for integration types which have one or two dimensions (z and/or y)
 * with limits of tau and 1, one more dimension with limits of 0 and 1 (xiprime),
 * and any number of other dimensions with limits plus or minus infinity
 */
class XiPIntegrationType : public IntegrationType {
public:
    XiPIntegrationType(size_t extra_dimensions) : IntegrationType(extra_dimensions) {}
protected:
    void fill_min(const Context& ctx, const size_t core_dimensions, double* min) const;
    void fill_max(const Context& ctx, const size_t core_dimensions, double* max) const;
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

/**
 * Base class for integration types which have one or two dimensions (z and/or y)
 * with limits of tau and 1, and any number of other dimensions in pairs with limits
 * Context::cutoff and Context::inf, and zero and 2pi, respectively.
 */
class RadialIntegrationType : public IntegrationType {
public:
    RadialIntegrationType(const size_t extra_dimensions) : IntegrationType(extra_dimensions) {assert(extra_dimensions % 2 == 0);}
    void fill_min(const Context& ctx, const size_t core_dimensions, double* min) const;
    void fill_max(const Context& ctx, const size_t core_dimensions, double* max) const;
};

class RadialPositionIntegrationType : public RadialIntegrationType {
public:
    RadialPositionIntegrationType(const size_t extra_dimensions) : RadialIntegrationType(extra_dimensions) {}
    double jacobian(const unsigned int ncoords, const double* coordinates, const IntegrationContext& ictx, const size_t core_dimensions) const;
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

class RadialMomentumIntegrationType : public RadialIntegrationType {
public:
    RadialMomentumIntegrationType(const size_t extra_dimensions) : RadialIntegrationType(extra_dimensions) {}
    double jacobian(const unsigned int ncoords, const double* coordinates, const IntegrationContext& ictx, const size_t core_dimensions) const;
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

class AngleIndependentPositionIntegrationType : public IntegrationType {
public:
    AngleIndependentPositionIntegrationType(const size_t extra_dimensions) : IntegrationType(extra_dimensions) {}
    void fill_min(const Context& ctx, const size_t core_dimensions, double *min) const;
    void fill_max(const Context& ctx, const size_t core_dimensions, double* max) const;
    double jacobian(const unsigned int ncoords, const double* coordinates, const IntegrationContext& ictx, const size_t core_dimensions) const;
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

class RescaledAngleIndependentPositionIntegrationType : public AngleIndependentPositionIntegrationType {
public:
    RescaledAngleIndependentPositionIntegrationType(const size_t extra_dimensions) : AngleIndependentPositionIntegrationType(extra_dimensions) {}
    double jacobian(const unsigned int ncoords, const double* coordinates, const IntegrationContext& ictx, const size_t core_dimensions) const;
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

class QLimitedMomentumIntegrationType : public RadialMomentumIntegrationType {
public:
    QLimitedMomentumIntegrationType(const size_t extra_dimensions) : RadialMomentumIntegrationType(extra_dimensions) {}
    void fill_max(const Context& ctx, const size_t core_dimensions, double* max) const;
    double jacobian(const unsigned int ncoords, const double* coordinates, const IntegrationContext& ictx, const size_t core_dimensions) const;
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

#endif
