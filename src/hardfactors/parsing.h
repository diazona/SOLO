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

/**
 * Attempts to parse a string as a hard factor group specification.
 * The specification passed should take the form
 *
 *     name:hfname.hfimpl,hfname.hfimpl,hfname.hfimpl,...
 *
 * where
 *
 * - `name` is the descriptive name of the hard factor group,  which can
 *   be used to refer to it on the command line of `oneloopcalc`. The name
 *   is optional; if one is not provided, the hard factor group's name will
 *   simply be the full text of the specification.
 * - `hfname.hfimpl` is a hard factor specification, as described in the
 *   documentation for ::parse_hardfactor.
 *
 * @return the group object, or `NULL` if the parsing failed to produce
 * a valid group
 */
const HardFactorGroup* parse_hardfactor_group(const std::string& spec);

/**
 * Splits a string representing a hard factor specification into a name
 * and implementation. The specification passed should take the form
 *
 *     name.implementation
 *
 * or just
 *
 *     name
 *
 * where
 *
 * - `name` is the descriptive name of the hard factor.
 * - `implementation` identifies which implementation of the hard factor
 *   is desired. Sometimes there are multiple expressions that can be
 *   used to calculate the same thing; the implementation code is a way
 *   to distinguish between them. For example, `m` or `momentum` means
 *   to use a momentum-space expression, `r` or `radial` for a radial
 *   position space expression, and `p` or `position` for a Cartesian
 *   position space expression. The implementation code is optional.
 *
 * This function makes no attempt to verify that either `name` or
 * `implementation` actually refers to a real hard factor. To do that,
 * use ::parse_hardfactor.
 *
 * @param[in] spec the hard factor specification
 * @param[out] name the name part of the specification
 * @param[out] implementation the implementation part of the specification,
 * or the empty string if no implementation part is provided
 */
void split_hardfactor(const std::string& spec, std::string& name, std::string& implementation);

/**
 * Attempts to parse a string as a hard factor specification.
 * The specification passed should take the form
 *
 *     name.implementation
 *
 * or just
 *
 *     name
 *
 * where
 *
 * - `name` is the descriptive name of the hard factor.
 * - `implementation` identifies which implementation of the hard factor
 *   is desired. Sometimes there are multiple expressions that can be
 *   used to calculate the same thing; the implementation code is a way
 *   to distinguish between them. For example, `m` or `momentum` means
 *   to use a momentum-space expression, `r` or `radial` for a radial
 *   position space expression, and `p` or `position` for a Cartesian
 *   position space expression. The implementation code is optional.
 *
 * This implementation of the function does not have backwards
 * compatibility with previous versions of SOLO.
 *
 * @param[in] spec the hard factor specification
 * @return the `HardFactor`, or `NULL` if no hard factor with the given
 * name and implementation could be found.
 */
const HardFactor* parse_hardfactor(const string& spec);

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
    InvalidHardFactorSpecException(const std::string& hfspec, const std::string& message) throw();
    InvalidHardFactorSpecException(const InvalidHardFactorSpecException& other) throw() : _message(other._message) {}
    ~InvalidHardFactorSpecException() throw() {}
    void operator=(const InvalidHardFactorSpecException& other);
    const char* what() const throw();
};
