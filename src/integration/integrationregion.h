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

#include <vector>
#include "integrationcontext.h"

/**
 * Superclass for ::IntegrationRegion components that handle integration
 * over `z` and `xi`.
 */
class CoreIntegrationRegion {
public:
    CoreIntegrationRegion();
    virtual ~CoreIntegrationRegion();

    /**
     * Writes the lower bounds of the integration region into the array
     * `min` using information from the `Context`.
     */
    virtual void fill_min(const Context& ctx, const bool xi_preintegrated_term, double* min) const = 0;
    /**
     * Writes the upper bounds of the integration region into the array
     * `max` using information from the `Context`.
     */
    virtual void fill_max(const Context& ctx, const bool xi_preintegrated_term, double* max) const = 0;
    /**
     * Calculates the jacobian determinant to convert from the variables being
     * integrated over to the variables involved in the calculation.
     *
     * This gets called after `update()` at each step of the calculation,
     * so the `IntegrationContext` has been updated with the current coordinates.
     */
    virtual double jacobian(const IntegrationContext& ictx, const bool xi_preintegrated_term) const = 0;
    /**
     * Updates the variables in the `IntegrationContext` based on the `values`,
     * which were passed in from the Monte Carlo integration routine.
     */
    virtual void update(IntegrationContext& ictx, const bool xi_preintegrated_term, const double* values) const = 0;
    /**
     * The number of variables being integrated over.
     */
    virtual size_t dimensions(const bool xi_preintegrated_term) const = 0;
    virtual bool operator<(const CoreIntegrationRegion& other) const;
};

/**
 * Superclass for ::IntegrationRegion components that handle integration
 * over positions and momenta.
 */
class AuxiliaryIntegrationRegion {
public:
    AuxiliaryIntegrationRegion(const size_t dimensions);
    virtual ~AuxiliaryIntegrationRegion();

    /**
     * Writes the lower bounds of the integration region into the array
     * `min` using information from the `Context`.
     */
    virtual void fill_min(const Context& ctx, double* min) const = 0;
    /**
     * Writes the upper bounds of the integration region into the array
     * `max` using information from the `Context`.
     */
    virtual void fill_max(const Context& ctx, double* max) const = 0;
    /**
     * Calculates the jacobian determinant to convert from the variables being
     * integrated over to the variables involved in the calculation.
     *
     * This gets called after `update()` at each step of the calculation,
     * so the `IntegrationContext` has been updated with the current coordinates.
     */
    virtual double jacobian(const IntegrationContext& ictx) const = 0;
    /**
     * Updates the variables in the `IntegrationContext` based on the `values`,
     * which were passed in from the Monte Carlo integration routine.
     */
    virtual void update(IntegrationContext& ictx, const double* values) const = 0;
    /**
     * The number of variables being integrated over.
     */
    virtual size_t dimensions() const;
    virtual bool operator<(const AuxiliaryIntegrationRegion& other) const;
protected:
    size_t m_dimensions;
};

/**
 * Represents the properties of an integration region: the integration
 * variable(s), the limits, and the jacobian.
 *
 * An instance of this class is responsible for setting the number
 * of variables being integrated over and effectively specifying which
 * variables those are. In particular, it needs to set bounds on each
 * variable, and also translate the Monte Carlo integration variables
 * (which must have constant bounds) to the kinematic variables (whose
 * bounds can depend on each other) by calculating the jacobian determinant
 * and setting the value of the kinematic variables using provided values
 * of the integration variables.
 *
 * An `IntegrationRegion` is made of components: one ::CoreIntegrationRegion,
 * which handles integration over `z` and `xi`, and zero or more
 * ::AuxiliaryIntegrationRegion instances, which handle integration over any
 * additional variables like positions or momenta.
 */
class IntegrationRegion {
public:
    IntegrationRegion(const CoreIntegrationRegion& core);
    IntegrationRegion(const CoreIntegrationRegion& core, const AuxiliaryIntegrationRegion& aux);
    IntegrationRegion(const CoreIntegrationRegion& core, const std::vector<const AuxiliaryIntegrationRegion*> auxregions);
    IntegrationRegion(const CoreIntegrationRegion& core, const size_t n_auxregions, const AuxiliaryIntegrationRegion** auxregions);
    virtual ~IntegrationRegion();

    /**
     * Writes the lower bounds of the integration region into the array
     * `min` using information from the `Context`.
     */
    virtual void fill_min(const Context& ctx, const bool xi_preintegrated_term, double* min) const;
    /**
     * Writes the upper bounds of the integration region into the array
     * `max` using information from the `Context`.
     */
    virtual void fill_max(const Context& ctx, const bool xi_preintegrated_term, double* max) const;
    /**
     * Calculates the jacobian determinant to convert from the variables being
     * integrated over to the variables involved in the calculation.
     *
     * This gets called after `update()` at each step of the calculation,
     * so the `IntegrationContext` has been updated with the current coordinates.
     */
    virtual double jacobian(const IntegrationContext& ictx, const bool xi_preintegrated_term) const;
    /**
     * Updates the variables in the `IntegrationContext` based on the `values`,
     * which were passed in from the Monte Carlo integration routine.
     */
    virtual void update(IntegrationContext& ictx, const bool xi_preintegrated_term, const double* values) const;
    /**
     * The number of variables being integrated over.
     */
    virtual size_t dimensions(const bool xi_preintegrated_term) const;

    virtual bool operator<(const IntegrationRegion& other) const;
    virtual IntegrationRegion operator*=(const AuxiliaryIntegrationRegion& other);

    const CoreIntegrationRegion& m_core_region;
    const std::vector<const AuxiliaryIntegrationRegion*> m_subregions;
};

/**
 * An integration subregion that integrates over z from tau to 1.
 */
class LOKinematicsIntegrationRegion : public CoreIntegrationRegion {
public:
    virtual void fill_min(const Context& ctx, const bool xi_preintegrated_term, double* min) const;
    virtual void fill_max(const Context& ctx, const bool xi_preintegrated_term, double* max) const;
    virtual double jacobian(const IntegrationContext& ictx, const bool xi_preintegrated_term) const;
    virtual void update(IntegrationContext& ictx, const bool xi_preintegrated_term, const double* values) const;
    virtual size_t dimensions(const bool xi_preintegrated_term) const;
};

/**
 * An integration subregion that integrates over z from tau to 1
 * and xi from tau/z to 1.
 */
class NLOKinematicsIntegrationRegion : public CoreIntegrationRegion {
public:
    /**
     * Constructs a new instance.
     *
     * @param lower_zero[in] whether to set the lower limit of the `xi`
     * integration to zero, or leave it at the default of `tau/z`
     */
    NLOKinematicsIntegrationRegion(const bool lower_zero = false);
    virtual void fill_min(const Context& ctx, const bool xi_preintegrated_term, double* min) const;
    virtual void fill_max(const Context& ctx, const bool xi_preintegrated_term, double* max) const;
    virtual double jacobian(const IntegrationContext& ictx, const bool xi_preintegrated_term) const;
    virtual void update(IntegrationContext& ictx, const bool xi_preintegrated_term, const double* values) const;
    virtual size_t dimensions(const bool xi_preintegrated_term) const;

private:
    const bool m_lower_zero;
};

/**
 * An integration subregion that integrates over z from tauhat to 1
 * and xi from tau/z to 1 - pT/(z*sqs)*exp(-y).
 */
class NLOClippedKinematicsIntegrationRegion : public CoreIntegrationRegion {
    // implementation note: this could perhaps be absorbed
    // into NLOKinematicsIntegrationRegion
    // with a few appropriate constructor parameters
public:
    virtual void fill_min(const Context& ctx, const bool xi_preintegrated_term, double* min) const;
    virtual void fill_max(const Context& ctx, const bool xi_preintegrated_term, double* max) const;
    virtual double jacobian(const IntegrationContext& ictx, const bool xi_preintegrated_term) const;
    virtual void update(IntegrationContext& ictx, const bool xi_preintegrated_term, const double* values) const;
    virtual size_t dimensions(const bool xi_preintegrated_term) const;
};

/**
 * A superclass for integration subregions parametrized in Cartesian
 * coordinates, ranging from -inf to inf.
 */
class CartesianIntegrationRegion : public AuxiliaryIntegrationRegion {
public:
    CartesianIntegrationRegion();
    virtual void fill_min(const Context& ctx, double* min) const;
    virtual void fill_max(const Context& ctx, double* max) const;
    virtual double jacobian(const IntegrationContext& ictx) const;
    virtual void update(IntegrationContext& ictx, const double* values) const;
protected:
    virtual double& x(IntegrationContext& ictx) const = 0;
    virtual double& y(IntegrationContext& ictx) const = 0;
    virtual double x(const IntegrationContext& ictx) const = 0;
    virtual double y(const IntegrationContext& ictx) const = 0;
};

/**
 * A superclass for integration subregions parametrized in polar
 * coordinates, with one ranging from 0 to inf and the other from
 * 0 to 2pi.
 *
 * The `m_center` field, initialized from the constructor, allows
 * for a shift in the origin of the polar coordinates.
 */
class PolarIntegrationRegion : public AuxiliaryIntegrationRegion {
public:
    PolarIntegrationRegion();
    virtual void fill_min(const Context& ctx, double* min) const;
    virtual void fill_max(const Context& ctx, double* max) const;
    virtual double jacobian(const IntegrationContext& ictx) const;
    virtual void update(IntegrationContext& ictx, const double* values) const;
protected:
    PolarIntegrationRegion(const size_t dimensions);
    virtual double& x(IntegrationContext& ictx) const = 0;
    virtual double& y(IntegrationContext& ictx) const = 0;
    virtual double x(const IntegrationContext& ictx) const = 0;
    virtual double y(const IntegrationContext& ictx) const = 0;
    virtual double r(const IntegrationContext& ictx) const;
    virtual void center(const IntegrationContext& ictx, double* center_point) const;
};

/**
 * A superclass for integration subregions parametrized by the
 * radial polar coordinate only. This is for integrands which are
 * independent of angle. The one coordinate ranges from 0 to inf.
 */
class RadialIntegrationRegion : public PolarIntegrationRegion {
public:
    RadialIntegrationRegion();
    virtual void fill_min(const Context& ctx, double* min) const;
    virtual void fill_max(const Context& ctx, double* max) const;
    virtual double jacobian(const IntegrationContext& ictx) const;
    virtual void update(IntegrationContext& ictx, const double* values) const;
};

/**
 * A superclass for integration subregions parametrized in polar
 * coordinates, with an explicit upper cutoff on the radial coordinate.
 */
class PolarLimitedIntegrationRegion : public PolarIntegrationRegion {
public:
    virtual void fill_max(const Context& ctx, double* max) const;
    virtual double jacobian(const IntegrationContext& ictx) const;
    virtual void update(IntegrationContext& ictx, const double* values) const;
protected:
    virtual double limit(const IntegrationContext& ictx) const = 0;
};

/**
 * This class implements the upper integration limit that emerges from
 * the exact kinematic constraint, which appears in e.g. arxiv:1405.6311
 * (the matching paper) and arxiv:1608.05293. The upper limit is
 *  q^2 < kT * (sqrt(s) * e^Y - kT) * (1 - xi) / xi
 * or equivalently,
 *  q^2 < (xp * s - kT^2) * (1 - xi) / xi
 */
class ExactKinematicIntegrationRegion : public PolarLimitedIntegrationRegion {
protected:
    virtual double limit(const IntegrationContext& ictx) const;
    virtual void center(const IntegrationContext& ictx, double* center_point) const;
};

/* Now we need versions of the cartesian, polar, and radial integration
 * regions for each of the three position coordinates and each
 * of the three momentum coordinates.
 */
template<class T>
class R1IntegrationRegion : public T {
protected:
    virtual double& x(IntegrationContext& ictx) const { return ictx.xx; }
    virtual double& y(IntegrationContext& ictx) const { return ictx.xy; }
    virtual double x(const IntegrationContext& ictx) const { return ictx.xx; }
    virtual double y(const IntegrationContext& ictx) const { return ictx.xy; }
};

template<class T>
class R2IntegrationRegion : public T {
protected:
    virtual double& x(IntegrationContext& ictx) const { return ictx.yx; }
    virtual double& y(IntegrationContext& ictx) const { return ictx.yy; }
    virtual double x(const IntegrationContext& ictx) const { return ictx.yx; }
    virtual double y(const IntegrationContext& ictx) const { return ictx.yy; }
};

template<class T>
class R3IntegrationRegion : public T {
protected:
    virtual double& x(IntegrationContext& ictx) const { return ictx.bx; }
    virtual double& y(IntegrationContext& ictx) const { return ictx.by; }
    virtual double x(const IntegrationContext& ictx) const { return ictx.bx; }
    virtual double y(const IntegrationContext& ictx) const { return ictx.by; }
};

template<class T>
class Q1IntegrationRegion : public T {
protected:
    virtual double& x(IntegrationContext& ictx) const { return ictx.q1x; }
    virtual double& y(IntegrationContext& ictx) const { return ictx.q1y; }
    virtual double x(const IntegrationContext& ictx) const { return ictx.q1x; }
    virtual double y(const IntegrationContext& ictx) const { return ictx.q1y; }
    virtual double r(const IntegrationContext& ictx) const { return sqrt(ictx.q12); }
};

template<class T>
class Q2IntegrationRegion : public T {
protected:
    virtual double& x(IntegrationContext& ictx) const { return ictx.q2x; }
    virtual double& y(IntegrationContext& ictx) const { return ictx.q2y; }
    virtual double x(const IntegrationContext& ictx) const { return ictx.q2x; }
    virtual double y(const IntegrationContext& ictx) const { return ictx.q2y; }
    virtual double r(const IntegrationContext& ictx) const { return sqrt(ictx.q22); }
};

template<class T>
class Q3IntegrationRegion : public T {
protected:
    virtual double& x(IntegrationContext& ictx) const { return ictx.q3x; }
    virtual double& y(IntegrationContext& ictx) const { return ictx.q3y; }
    virtual double x(const IntegrationContext& ictx) const { return ictx.q3x; }
    virtual double y(const IntegrationContext& ictx) const { return ictx.q3y; }
    virtual double r(const IntegrationContext& ictx) const { return sqrt(ictx.q32); }
};

/* This creates specific subclasses of e.g. CartesianIntegrationRegion
 * for each of R1, R2, R3, Q1, Q2, and Q3. The names will be e.g.
 *  R1CartesianIntegrationRegion
 *  R2CartesianIntegrationRegion
 *  R3CartesianIntegrationRegion
 *  Q1CartesianIntegrationRegion
 *  Q2CartesianIntegrationRegion
 *  Q3CartesianIntegrationRegion
 */
#define MAKE_CLASSES(T) \
typedef R1IntegrationRegion<T> R1 ## T;\
typedef R2IntegrationRegion<T> R2 ## T;\
typedef R3IntegrationRegion<T> R3 ## T;\
typedef Q1IntegrationRegion<T> Q1 ## T;\
typedef Q2IntegrationRegion<T> Q2 ## T;\
typedef Q3IntegrationRegion<T> Q3 ## T;

MAKE_CLASSES(CartesianIntegrationRegion)
MAKE_CLASSES(PolarIntegrationRegion)
MAKE_CLASSES(RadialIntegrationRegion)
MAKE_CLASSES(ExactKinematicIntegrationRegion)

#undef MAKE_CLASSES

#endif
