/*
 * Part of SOLO
 * 
 * Copyright 2014 David Zaslavsky
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
#include "parsing.h"
#include "../utils/utils.h"

void split_hardfactor(const string& spec, string& name, string& implementation) {
    vector<string> splitspec;
    splitspec = split(spec, ".", 2);
    assert(splitspec.size() == 1 || splitspec.size() == 2);
    name = splitspec[0];
    if (splitspec.size() == 2) {
        implementation = splitspec[1];
    }
    else {
        implementation = "";
    }
}

const HardFactor* parse_hardfactor(const string& spec) {
    // the hard factor takes the form <name>.<implementation>
    string name;
    string implementation;
    split_hardfactor(spec, name, implementation);
    const HardFactor* hardfactor;
    if (!implementation.empty()) {
        if (implementation == "momentum" || implementation == "m") {
            return momentum::registry::get_instance()->get_hard_factor(name);
        }
        else if (implementation == "radial" || implementation == "r") {
            return radial::registry::get_instance()->get_hard_factor(name);
        }
        else if (implementation == "position" || implementation == "p") {
            return position::registry::get_instance()->get_hard_factor(name);
        }
        else {
            return NULL;
        }
    }
    else {
        // by default try momentum first, then radial, then position
        hardfactor = momentum::registry::get_instance()->get_hard_factor(name);
        if (hardfactor != NULL) {
            return hardfactor;
        }
        hardfactor = radial::registry::get_instance()->get_hard_factor(name);
        if (hardfactor != NULL) {
            return hardfactor;
        }
        hardfactor = position::registry::get_instance()->get_hard_factor(name);
        return hardfactor; // whether NULL or not
    }
}

const HardFactorGroup* parse_hardfactor_group(const string& spec) {
    string specname(spec);
    string specbody(spec);

    vector<string> namesplitspec = split(spec, ":", 2);
    assert(namesplitspec.size() == 1 || namesplitspec.size() == 2);
    if (namesplitspec.size() > 1) {
        if (namesplitspec[0].size() == 0) {
            // this means there was simply a leading colon; ignore it
            specname = namesplitspec[1];
        }
        else {
            specname = namesplitspec[0];
        }
        specbody = namesplitspec[1];
    }

    vector<string> splitspec;
    splitspec = split(specbody, ", ");

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

InvalidHardFactorSpecException::InvalidHardFactorSpecException(const string& hfspec, const string& message) throw() : hfspec(hfspec) {
    std::ostringstream s;
    s << message << "; spec " << hfspec;
    _message = s.str();
}

void InvalidHardFactorSpecException::operator=(const InvalidHardFactorSpecException& other) {
    _message = other._message;
}

const char* InvalidHardFactorSpecException::what() const throw() {
    return _message.c_str();
}
