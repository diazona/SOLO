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
 * Represents a particular integration region: a combination of the
 * number of dimensions being integrated over and the bounds in each
 * of those dimensions.
 */
class IntegrationType {
public:
    /** The number of variables being integrated over excluding z and y */
    const size_t extra_dimensions;
    IntegrationType(const size_t extra_dimensions) : extra_dimensions(extra_dimensions) {}
    /**
     * Writes the lower bounds of the integration region into the array
     * `min` using information from the `IntegrationContext`.
     */
    virtual void fill_min(IntegrationContext& ictx, const size_t core_dimensions, double* min) const = 0;
    /**
     * Writes the upper bounds of the integration region into the array
     * `max` using information from the `IntegrationContext`.
     */
    virtual void fill_max(IntegrationContext& ictx, const size_t core_dimensions, double* max) const = 0;
    /**
     * Updates the variables in the `IntegrationContext` based on the `values`,
     * which were passed in from the Monte Carlo integration routine.
     */
    virtual void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const = 0;
};

/**
 * Base class for integration types which have one or two dimensions (z and/or y)
 * with limits of tau and 1, and any number of other dimensions with limits
 * plus or minus infinity
 */
class PlainIntegrationType : public IntegrationType {
protected:
    PlainIntegrationType(const size_t extra_dimensions) : IntegrationType(extra_dimensions) {}
    void fill_min(IntegrationContext& ictx, const size_t core_dimensions, double* min) const;
    void fill_max(IntegrationContext& ictx, const size_t core_dimensions, double* max) const;
};

/**
 * An integration type for integration over rx, ry, z and/or y
 */
class DipoleIntegrationType : public PlainIntegrationType {
public:
    static const DipoleIntegrationType* get_instance() {
        static DipoleIntegrationType instance;
        return &instance;
    }
private:
    DipoleIntegrationType() : PlainIntegrationType(2) {}
    DipoleIntegrationType(DipoleIntegrationType const&);
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

/**
 * An integration type for integration over sx, sy, tx, ty, z and/or y
 */
class QuadrupoleIntegrationType : public PlainIntegrationType {
public:
    static const QuadrupoleIntegrationType* get_instance() {
        static QuadrupoleIntegrationType instance;
        return &instance;
    }
private:
    QuadrupoleIntegrationType() : PlainIntegrationType(4) {}
    QuadrupoleIntegrationType(QuadrupoleIntegrationType const&);
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

/**
 * An integration type for integration over q1x, q1y, z and/or y
 */
class Momentum1IntegrationType : public PlainIntegrationType {
public:
    static const Momentum1IntegrationType* get_instance() {
        static Momentum1IntegrationType instance;
        return &instance;
    }
private:
    Momentum1IntegrationType() : PlainIntegrationType(2) {}
    Momentum1IntegrationType(Momentum1IntegrationType const&);
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

/**
 * An integration type for integration over q1x, q1y, q2x, q2y, z and/or y
 */
class Momentum2IntegrationType : public PlainIntegrationType {
public:
    static const Momentum2IntegrationType* get_instance() {
        static Momentum2IntegrationType instance;
        return &instance;
    }
private:
    Momentum2IntegrationType() : PlainIntegrationType(4) {}
    Momentum2IntegrationType(Momentum2IntegrationType const&);
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

/**
 * An integration type for integration over q1x, q1y, q2x, q2y, q3x, q3y, z and/or y
 */
class Momentum3IntegrationType : public PlainIntegrationType {
public:
    static const Momentum3IntegrationType* get_instance() {
        static Momentum3IntegrationType instance;
        return &instance;
    }
private:
    Momentum3IntegrationType() : PlainIntegrationType(6) {}
    Momentum3IntegrationType(Momentum3IntegrationType const&);
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

/**
 * Base class for integration types which have one or two dimensions (z and/or y)
 * with limits of tau and 1, one more dimension with limits of 0 and 1 (xiprime),
 * and any number of other dimensions with limits plus or minus infinity
 */
class XiPIntegrationType : public IntegrationType {
protected:
    XiPIntegrationType(size_t extra_dimensions) : IntegrationType(extra_dimensions) {}
    void fill_min(IntegrationContext& ictx, const size_t core_dimensions, double* min) const;
    void fill_max(IntegrationContext& ictx, const size_t core_dimensions, double* max) const;
};

/**
 * An integration type for integration over q1x, q1y, xiprime, z and/or y
 */
class Momentum1XiPIntegrationType : public XiPIntegrationType {
public:
    static const Momentum1XiPIntegrationType* get_instance() {
        static Momentum1XiPIntegrationType instance;
        return &instance;
    }
private:
    Momentum1XiPIntegrationType() : XiPIntegrationType(3) {}
    Momentum1XiPIntegrationType(Momentum1XiPIntegrationType const&);
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

/**
 * An integration type for integration over q1x, q1y, q2x, q2y, xiprime, z and/or y
 */
class Momentum2XiPIntegrationType : public XiPIntegrationType {
public:
    static const Momentum2XiPIntegrationType* get_instance() {
        static Momentum2XiPIntegrationType instance;
        return &instance;
    }
private:
    Momentum2XiPIntegrationType() : XiPIntegrationType(5) {}
    Momentum2XiPIntegrationType(Momentum2XiPIntegrationType const&);
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

extern const double inf;

#endif
