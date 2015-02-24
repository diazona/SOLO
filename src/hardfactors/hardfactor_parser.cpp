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
#include "hardfactor_parser.h"

HardFactorRegistry* parse_one_hardfactor_spec(const string& spec, string& name, const bool check_for_existing) {
    char hf_type = '\000'; // dummy value, as a default
    
    // If the hard factor specification takes the form [letter]:[stuff],
    // consider the first letter to indicate the type, and remove it and
    // the colon
    if (spec[1] == '.' || spec[1] == ':') {
        hf_type = spec[0];
        if (is_registry(hf_type)) {
            name = spec.substr(2);
        }
        else {
            hf_type = '\000';
        }
    }
    // Pass the remaining name (e.g. "h02qq") to the hard factor registry
    // to get the actual hard factor object
    switch (hf_type) {
        case 'm':
            return momentum::registry::get_instance();
        case 'r':
            return radial::registry::get_instance();
        case 'p':
            return position::registry::get_instance();
        default:
            assert(!is_registry(hf_type));
            if (check_for_existing) {
                // by default try momentum first, then radial, then position
                const HardFactor* hf = momentum::registry::get_instance()->get_hard_factor(spec);
                if (hf != NULL) {
                    name = spec;
                    return momentum::registry::get_instance();
                }
                hf = radial::registry::get_instance()->get_hard_factor(spec);
                if (hf != NULL) {
                    name = spec;
                    return radial::registry::get_instance();
                }
                hf = position::registry::get_instance()->get_hard_factor(spec);
                if (hf != NULL) {
                    name = spec;
                    return position::registry::get_instance();
                }
            }
    }
    return NULL;
}


const HardFactor* parse_hardfactor(const string& spec) {
    string name;

    const HardFactorRegistry* registry = parse_one_hardfactor_spec(spec, name, true);
    return registry == NULL ? NULL : registry->get_hard_factor(name);
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
