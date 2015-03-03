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

#include "hardfactor.h"
#include "../utils/utils.h"

namespace position {
    const PositionIntegrationType dipole(2);
    const PositionIntegrationType quadrupole(4);
}

namespace radial {
    const AngleIndependentPositionIntegrationType dipole(1);
    const AngleIndependentPositionIntegrationType quadrupole(2);
    const RescaledAngleIndependentPositionIntegrationType rescaled_dipole(1);
    const RescaledAngleIndependentPositionIntegrationType rescaled_quadrupole(2);
}

namespace momentum {
    const NoIntegrationType none;
    const MomentumIntegrationType momentum1(2);
    const MomentumIntegrationType momentum2(4);
    const MomentumIntegrationType momentum3(6);
    const RadialMomentumIntegrationType radialmomentum1(2);
    const RadialMomentumIntegrationType radialmomentum2(4);
    const RadialMomentumIntegrationType radialmomentum3(6);
    const XiPIntegrationType momentumxip1(3);
    const XiPIntegrationType momentumxip2(5);
    const QLimitedMomentumIntegrationType qlim(2);
}

using std::string;
using std::vector;

// for weird technical reasons, this needs to be in a .cpp file, not a header
HardFactor::~HardFactor() {}

const HardFactor::HardFactorOrder HardFactor::get_order() const {
    /* This default implementation relies on a particular convention
     * for get_name(), but can be overridden for hard factors where
     * that convention doesn't apply.
     */
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
    else {
        assert(false);
    }
};

const size_t HardFactorTerm::get_term_count() const {
    return 1;
}

const HardFactorTerm* const* HardFactorTerm::get_terms() const {
    return &p_this;
}

void HardFactorRegistry::add_hard_factor(const HardFactor* hf, bool manage) {
    add_hard_factor(hf->get_name(), hf, manage);
}

void HardFactorRegistry::add_hard_factor(const char* key, const HardFactor* hf, bool manage) {
    using namespace std;
    string keystring(key);
    assert(!keystring.empty());
    transform(keystring.begin(), keystring.end(), keystring.begin(), ::tolower);
    // note that even if another HardFactor is later added with the
    // same name, it doesn't cause a memory leak
    hardfactors[keystring] = hf;
    if (manage) {
        hardfactors_to_delete.push_back(hf);
    }
}

const HardFactor* HardFactorRegistry::get_hard_factor(const string& key) const {
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

void HardFactorRegistry::add_hard_factor_group(const HardFactorGroup* hfg, bool manage) {
    add_hard_factor_group(hfg->label.c_str(), hfg, manage);
}

void HardFactorRegistry::add_hard_factor_group(const char* key, const HardFactorGroup* hfg, bool manage) {
    using namespace std;
    string keystring(key);
    assert(!keystring.empty());
    transform(keystring.begin(), keystring.end(), keystring.begin(), ::tolower);
    // note that even if another HardFactorGroup is later added with the
    // same name, it doesn't cause a memory leak
    hardfactor_groups[keystring] = hfg;
    if (manage) {
        hardfactor_groups_to_delete.push_back(hfg);
    }
}

const HardFactorGroup* HardFactorRegistry::get_hard_factor_group(const string& key) const {
    using namespace std;
    string keystring(key);
    transform(keystring.begin(), keystring.end(), keystring.begin(), ::tolower);
    map<const string, const HardFactorGroup*>::const_iterator it = hardfactor_groups.find(key);
    if (it == hardfactor_groups.end()) {
        return NULL;
    }
    else {
        return it->second;
    }
}

HardFactorRegistry::~HardFactorRegistry() {
    hardfactors.clear();
    hardfactor_groups.clear();
    for (std::list<const HardFactor*>::iterator it = hardfactors_to_delete.begin(); it != hardfactors_to_delete.end(); it++) {
        delete (*it);
    }
    hardfactors_to_delete.clear();
    for (std::list<const HardFactorGroup*>::iterator it = hardfactor_groups_to_delete.begin(); it != hardfactor_groups_to_delete.end(); it++) {
        delete (*it);
    }
    hardfactor_groups_to_delete.clear();
}

KinematicSchemeMismatchException::KinematicSchemeMismatchException(const HardFactor& hf) throw() {
    _message = "Mixed-order hard factor ";
    _message += hf.get_name();
    _message += " cannot be integrated in exact kinematics";
}

void KinematicSchemeMismatchException::operator=(const KinematicSchemeMismatchException& other) {
    _message = other._message;
}

const char* KinematicSchemeMismatchException::what() const throw() {
    return _message.c_str();
}
