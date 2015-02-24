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
#include "hardfactor_parser.h"
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

/**
 * This defines the hard factor group that is used when "lo" is given
 * on the command line
 */
static const char* default_lo_spec = "m.h02qq,m.h02gg";
/**
 * This defines the hard factor group that is used when "nlo" is given
 * on the command line and exact_kinematics is not set, or when "nlo.std"
 * is given on the command line
 */
static const char* standard_nlo_spec = "r.h12qq,m.h14qq,r.h12gg,m.h12qqbar,m.h16gg,r.h112gq,r.h122gq,m.h14gq,r.h112qg,r.h122qg,m.h14qg";
/**
 * This defines the hard factor group that is used when "nlo" is given
 * on the command line and exact_kinematics is set, or when "nlo.hipt"
 * is given on the command line
 */
static const char* highpt_nlo_spec = "m.h1qqexact,m.h1ggexact,r.h112gq,r.h122gq,m.h14gq,r.h112qg,r.h122qg,m.h14qg";

const HardFactorGroup* HardFactorGroup::parse(const string& spec, bool exact_kinematics) {
    vector<string> splitspec;
    string specname(spec);
    // Split the specification string on commas to get individual hard factor names
    if (spec == "lo") {
        splitspec = split(default_lo_spec, ", ");
    }
    else if (spec == "nlo") {
        splitspec = split(exact_kinematics ? highpt_nlo_spec : standard_nlo_spec, ", ");
    }
    else if (spec == "nlo.hipt") {
        splitspec = split(highpt_nlo_spec, ", ");
    }
    else if (spec == "nlo.std") {
        splitspec = split(standard_nlo_spec, ", ");
    }
    else {
        string specbody(spec);
        vector<string> namesplitspec = split(spec, ":", 2);
        assert(namesplitspec.size() == 1 || namesplitspec.size() == 2);
        if (namesplitspec.size() > 1) {
            if (namesplitspec[0].size() == 0) {
                // this means there was simply a leading colon; ignore it
                specname = namesplitspec[1];
                specbody = namesplitspec[1];
            }
            else {
                /* At this point we have established that the hf spec takes the form
                 * <string1>:<string2> where <string1> is of nonzero length and contains
                 * no colons. This next bit of logic attempts to intelligently determine
                 * whether <string1> is meant to be a hard factor type specification
                 * (i.e. one of 'm', 'r', or 'p') or a label for the hard factor spec.
                 */
                if (!is_registry(namesplitspec[0])) {
                    // <string1> is something other than 'm', 'r', or 'p', so is clearly a label
                    specname = namesplitspec[0];
                    specbody = namesplitspec[1];
                }
                else if (
                    (namesplitspec[1][1] == ':' || namesplitspec[1][1] == '.')
                    && is_registry(namesplitspec[0])
                ) {
                    // <string2> starts with [mrp][:.] (in regex syntax) so
                    // includes a type letter, thus <string1> must be a label
                    specname = namesplitspec[0];
                    specbody = namesplitspec[1];
                }
                // else specname and specbody should both just be spec, but
                // they already are, so do nothing
            }
        }
        splitspec = split(specbody, ", ");
    }
    HardFactorList hfobjs;
    // Iterate over the individual hard factor names
    for (vector<string>::iterator it = splitspec.begin(); it != splitspec.end(); it++) {
        const HardFactor* hf = parse_hardfactor(*it);
        if (hf == NULL) {
            // the string failed to parse
            throw InvalidHardFactorSpecException(spec, "No such hard factor");
        }
        else {
            hfobjs.push_back(hf);
        }
    }
    if (hfobjs.size() == 0) {
        throw InvalidHardFactorSpecException(spec, "No valid hard factors in specification");
    }
    else {
        HardFactorList* p_hfobjs = new HardFactorList(hfobjs);
        // parsing seems to have succeeded, so go ahead and create the object
        return new HardFactorGroup(specname, p_hfobjs, splitspec);
    }
}

HardFactorRegistry::~HardFactorRegistry() {
    hardfactors.clear();
    for (std::list<const HardFactor*>::iterator it = hardfactors_to_delete.begin(); it != hardfactors_to_delete.end(); it++) {
        delete (*it);
    }
    hardfactors_to_delete.clear();
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
