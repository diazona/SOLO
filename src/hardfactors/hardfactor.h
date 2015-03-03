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

#ifndef _HARD_FACTOR_H_
#define _HARD_FACTOR_H_

#include <algorithm>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "../typedefs.h"
#include "../integration/integrationcontext.h"
#include "../integration/integrationtype.h"

class HardFactorTerm;

/**
 * Something that can be integrated using the program.
 * 
 * A `HardFactor` is a named collection of one or more ::HardFactorTerm
 * instances. When a `HardFactor` is integrated, the program
 * queries it for its `HardFactorTerm`s and sorts those terms out by
 * their `IntegrationType`. It then iterates through the `IntegrationType`s
 * and for each type, integrates the terms. The integrand for a given
 * type is just the sum of values of all the terms of that type.
 * 
 * In practice, what is usually integrated is a hard factor group,
 * which is a list of multiple `HardFactor` objects. The procedure is
 * the same; the program queries all the `HardFactor`s for their terms
 * and puts them all together into one big map of `IntegrationType`s to sets
 * of `HardFactorTerm`s. The difference is that, when the `--separate` option
 * is passed to the main program, the `HardFactor`s in a hard factor group
 * will be integrated separately. But the `HardFactorTerm`s within a single
 * `HardFactor` are always kept together.
 */
class HardFactor {
public:
    typedef enum {LO, NLO, MIXED} HardFactorOrder;
    virtual ~HardFactor();
    /** The human-readable name of the hard factor. */
    virtual const char* get_name() const = 0;
    /** The number of HardFactorTerm objects this hard factor has. */
    virtual const size_t get_term_count() const = 0;
    /** A pointer to the list of HardFactorTerm objects. */
    virtual const HardFactorTerm* const* get_terms() const = 0;
    /**
     * The order of the term (LO, NLO, mixed). This sets which definitions
     * of kinematic variables (IntegrationContext::xp, IntegrationContext::xg)
     * are used when evaluating the `HardFactor`.
     */
    virtual const HardFactorOrder get_order() const;
};

/**
 * The base of the classes that actually implement the formulas.
 * 
 * A HardFactorTerm can be queried for three functions: `Fs`, `Fn`, `Fd`.
 * Each returns a real and an imaginary component. These three functions
 * are used in particular combinations to compute the "1D" and "2D" integrands.
 * (The mathematical details are explained elsewhere.)
 * 
 * A `HardFactorTerm` is also a `HardFactor` which contains just one term: itself.
 */
class HardFactorTerm : public HardFactor {
public:
    HardFactorTerm() : p_this(this) {};
    /** The IntegrationType needed to integrate this term. */
    virtual const IntegrationType* get_type() const = 0;
    /** The plus-regulated ("singular") part of the term. */
    virtual void Fs(const IntegrationContext* ictx, double* real, double* imag) const { *real = 0; *imag = 0; }
    /** The normal part of the term. */
    virtual void Fn(const IntegrationContext* ictx, double* real, double* imag) const { *real = 0; *imag = 0; }
    /** The delta-function part of the term. */
    virtual void Fd(const IntegrationContext* ictx, double* real, double* imag) const { *real = 0; *imag = 0; }
    /**
     * @return 1
     */
    const size_t get_term_count() const;
    /**
     * @return `this`
     */
    const HardFactorTerm* const* get_terms() const;
private:
    // necessary because `this` isn't an actual pointer that is stored anywhere
    const HardFactorTerm* p_this;
};

/**
 * A group of `HardFactor`s which can be calculated together or separately.
 */
class HardFactorGroup {
public:
    HardFactorGroup(const std::string& label, const HardFactorList* objects, const std::vector<std::string>& specifications) :
     label(label), objects(objects), specifications(specifications) {}
    const std::string label;
    const HardFactorList* objects;
    const std::vector<std::string> specifications;
};


/**
 * A class that holds a map of strings to HardFactor objects,
 * and can return the HardFactor object corresponding to a given string key.
 */
class HardFactorRegistry {
public:
    /**
     * Adds a HardFactor pointer to the registry. It can later be
     * retrieved by its name.
     * 
     * @param manage true if the HardFactor object should be deleted
     * by the registry
     */
    void add_hard_factor(const HardFactor* hf, bool manage=false);
    /**
     * Adds a HardFactor pointer to the registry under a custom key.
     * It can later be retrieved by the provided key.
     * 
     * @param manage true if the HardFactor object should be deleted
     * by the registry
     */
    void add_hard_factor(const char* key, const HardFactor* hf, bool manage=false);
    /**
     * Returns the hard factor corresponding to the given key, or NULL if none
     */
    const HardFactor* get_hard_factor(const std::string& key) const;
protected:
    ~HardFactorRegistry();
private:
    std::map<const std::string, const HardFactor*> hardfactors;
    std::list<const HardFactor*> hardfactors_to_delete;
};

/**
 * Exception to throw when trying to integrate a mixed LO+NLO hard factor
 * in exact kinematics. This is not possible because LO and NLO terms use
 * different integration limits in exact kinematics.
 */
class KinematicSchemeMismatchException : public std::exception {
private:
    string _message;
public:
    KinematicSchemeMismatchException(const HardFactor& hf) throw();
    KinematicSchemeMismatchException(const KinematicSchemeMismatchException& other) throw() : _message(other._message) {}
    ~KinematicSchemeMismatchException() throw() {}
    void operator=(const KinematicSchemeMismatchException& other);
    const char* what() const throw();
};

namespace position {

extern const PositionIntegrationType dipole;
extern const PositionIntegrationType quadrupole;

class registry : public HardFactorRegistry {
public:
    static registry* get_instance() {
        static registry instance;
        return &instance;
    }
private:
    registry() {}
};

}

namespace radial {

extern const AngleIndependentPositionIntegrationType dipole;
extern const AngleIndependentPositionIntegrationType quadrupole;
extern const RescaledAngleIndependentPositionIntegrationType rescaled_dipole;
extern const RescaledAngleIndependentPositionIntegrationType rescaled_quadrupole;

class registry : public HardFactorRegistry {
public:
    static registry* get_instance() {
        static registry instance;
        return &instance;
    }
private:
    registry() {}
};

}

namespace momentum {

extern const NoIntegrationType none;
extern const MomentumIntegrationType momentum1;
extern const MomentumIntegrationType momentum2;
extern const MomentumIntegrationType momentum3;
extern const RadialMomentumIntegrationType radialmomentum1;
extern const RadialMomentumIntegrationType radialmomentum2;
extern const RadialMomentumIntegrationType radialmomentum3;
extern const XiPIntegrationType momentumxip1;
extern const XiPIntegrationType momentumxip2;
extern const QLimitedMomentumIntegrationType qlim;

class registry : public HardFactorRegistry {
public:
    static registry* get_instance() {
        static registry instance;
        return &instance;
    }
private:
    registry() {}
};

}


#endif // _HARD_FACTOR_H_
