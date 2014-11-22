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

#pragma once

#include <stdexcept>
#include "hardfactor.h"

const HardFactor* parse_hardfactor(const string& spec);

/**
 * Splits a specification for a single hard factor into a registry and a name.
 * 
 * @param[in] spec the specification, `m.h02qq` or similar
 * @param[out] name reference to a variable that will be set to the name,
 * `h02qq` in this case. However, if this function returns `NULL`, the value
 * of `name` will not be changed.
 * @param[in] check_for_existing if the spec doesn't identify a registry, then
 * if this variable is `true`, the function will check all the registries for
 * a hard factor with that name and return the first one it finds. If this
 * variable is `false` and the spec doesn't identify a registry, the function
 * will simply return `NULL`. This is useful to set to `true` when searching
 * for an existing hard factor, and `false` when parsing a new hard factor.
 * @return a pointer to the registry identified by the spec, or the registry
 * in which a hard factor with the given name was found, if that check was
 * performed.
 */
HardFactorRegistry* parse_one_hardfactor_spec(const string& spec, string& name, const bool check_for_existing);

inline bool is_registry(const char implementation_code) {
    switch (implementation_code) {
        case 'm':
        case 'r':
        case 'p':
            return true;
        default:
            return false;
    }
}

inline bool is_registry(const string& implementation_code) {
    return implementation_code.length() == 1 ? is_registry(implementation_code[0]) : false;
}

/**
 * An exception to be thrown when a hard factor specification fails to be
 * parsed for some reason. The reason is indicated by the message.
 */
class InvalidHardFactorSpecException : public std::exception {
private:
    std::string _message;
public:
    std::string hfspec;
    /**
     * Constructs an instance of the exception.
     *
     * @param hfspec the hard factor specification or part of a specification
     *  that caused the error
     * @param message a descriptive human-readable message indicating why hfspec
     *  could not be parsed
     */
    InvalidHardFactorSpecException(const std::string& hfspec, const std::string& message) throw() : hfspec(hfspec) {
        std::ostringstream s;
        s << message << " in hard factor " << hfspec;
        _message = s.str();
    }
    InvalidHardFactorSpecException(const InvalidHardFactorSpecException& other) throw() : _message(other._message) {}
    ~InvalidHardFactorSpecException() throw() {}
    void operator=(const InvalidHardFactorSpecException& other) {
        _message = other._message;
    }
    const char* what() const throw() {
        return _message.c_str();
    }
};
