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
#include <string>
#include "integrationcontext.h"
#include "integrationtype.h"

class HardFactorTerm;

/**
 * Something that can be integrated using the program.
 * 
 * The actual implementation of the formulas is left for subclasses
 * of HardFactorTerm. A HardFactor can be either a single HardFactorTerm,
 * or a group of them. When a HardFactor is integrated, the program
 * queries it for its HardFactorTerms and sorts those terms out by
 * their IntegrationType. It then iterates through the IntegrationTypes
 * and for each type, integrates the terms. The integrand for a given
 * type is just the sum of values of all the terms of that type.
 * 
 * In practice, what is usually integrated is a hard factor group,
 * which is a list of multiple HardFactor objects. The procedure is
 * the same, it just puts all the HardFactorTerms together into one
 * big map.
 */
class HardFactor {
public:
    typedef enum {LO, NLO, MIXED} HardFactorOrder;
    /** An identifying name for the hard factor. */
    virtual const char* get_name() const = 0;
    /** The number of HardFactorTerm objects this hard factor has. */
    virtual const size_t get_term_count() const = 0;
    /** A pointer to the list of HardFactorTerm objects. */
    virtual const HardFactorTerm* const* get_terms() const = 0;
    /** The order of the term (LO, NLO, mixed) */
    virtual HardFactorOrder get_order() const {
        // relies on a particular convention for get_name()
        // but can be overridden for hard factors where that convention doesn't apply
        std::string name = get_name();
        if (name.compare(0, 3, "H01") == 0) {
            return MIXED;
        }
        else if (name.compare(0, 2, "H0") == 0) {
            return LO;
        }
        else if (name.compare(0, 2, "H1") == 0) {
            return NLO;
        }
    };
};

/**
 * The base of the classes that actually implement the formulas.
 * 
 * A HardFactorTerm can be queried for three functions: Fs, Fn, Fd.
 * Each returns a real and an imaginary component. These three functions
 * are used in particular combinations to compute the "1D" and "2D" integrands.
 * (The mathematical details are explained elsewhere.)
 * 
 * A HardFactorTerm is also a HardFactor which contains just one term: itself.
 */
class HardFactorTerm : public HardFactor {
public:
    HardFactorTerm() : p_this(this) {};
    /** Returns the IntegrationType needed to integrate this term. */
    virtual const IntegrationType* get_type() const = 0;
    virtual void Fs(const IntegrationContext* ictx, double* real, double* imag) const { *real = 0; *imag = 0; }
    virtual void Fn(const IntegrationContext* ictx, double* real, double* imag) const { *real = 0; *imag = 0; }
    virtual void Fd(const IntegrationContext* ictx, double* real, double* imag) const { *real = 0; *imag = 0; }
    const size_t get_term_count() const {
        return 1;
    }
    const HardFactorTerm* const* get_terms() const {
        return &p_this;
    }
private:
    const HardFactorTerm* p_this; // wtf why do I have to do this
};

/**
 * A singleton class that holds a map of strings to HardFactor objects,
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
    void add_hard_factor(const HardFactor* hf, bool manage=false) {
        add_hard_factor(hf->get_name(), hf, manage);
    }
    /**
     * Adds a HardFactor pointer to the registry under a custom key.
     * It can later be retrieved by the provided key.
     * 
     * @param manage true if the HardFactor object should be deleted
     * by the registry
     */
    void add_hard_factor(const char* key, const HardFactor* hf, bool manage=false) {
        using namespace std;
        string keystring(key);
        transform(keystring.begin(), keystring.end(), keystring.begin(), ::tolower);
        hardfactors[keystring] = hf;
        if (manage) {
            hardfactors_to_delete.push_back(hf);
        }
    }
    /**
     * Returns the hard factor corresponding to the given key, or NULL if none
     */
    const HardFactor* get_hard_factor(std::string& key) const {
        using namespace std;
        string keystring(key);
        transform(keystring.begin(), keystring.end(), keystring.begin(), ::tolower);
        map<const string, const HardFactor*>::const_iterator it = hardfactors.find(key);
        if (it == hardfactors.end()) {
            return NULL;
        }
        else {
            return it->second;
        }
    }
protected:
    ~HardFactorRegistry() {
        hardfactors.clear();
        for (std::list<const HardFactor*>::iterator it = hardfactors_to_delete.begin(); it != hardfactors_to_delete.end(); it++) {
            delete (*it);
        }
        hardfactors_to_delete.clear();
    }
private:
    std::map<const std::string, const HardFactor*> hardfactors;
    std::list<const HardFactor*> hardfactors_to_delete;
};

#endif // _HARD_FACTOR_H_
