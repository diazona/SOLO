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
#include <sstream>
#include <string>
#include <vector>
#include "../typedefs.h"

/**
 * A class to parse a hard factor specification.
 */
class ParsedHardFactorGroup {
private:
    ParsedHardFactorGroup(const std::string& label, const HardFactorList* objects, const std::vector<std::string>& specifications) :
     label(label), objects(objects), specifications(specifications) {}
public:
    /**
     * Attempts to parse a string as a hard factor group specification.
     * 
     * The specification passed should take the form
     * 
     *     name:x.hfname,x.hfname,x.hfname,...
     * 
     * where
     * 
     * - `name` is the descriptive name of the hard factor group,  which can
     *   be used to refer to it on the command line of `oneloopcalc`. The name
     *   is optional; if one is not provided, the hard factor group's name will
     *   simply be the full text of the specification.
     * - `x` is an implementation code, which identifies which ::HardFactorRegistry
     *   will be searched for the hard factor. For example, `m` identifies the
     *   registry for momentum-space hard factors, momentum::registry.
     *   Similarly `r` is for radial::registry and `p` for position::registry.
     * - `hfname` is the descriptive name of the hard factor itself. The registry
     *   identified by the implementation code will be queried for the hard factor
     *   by passing `hfname` to the HardFactorRegistry::get_hard_factor() method.
     * 
     * For backwards compatibility, the period in `x.hfname` can be a colon
     * instead.
     *
     * Returns `NULL` if the parsing failed to produce a valid group.
     */
    static const ParsedHardFactorGroup* parse(const std::string& spec, bool exact_kinematics);

    const std::string label;
    const HardFactorList* objects;
    const std::vector<std::string> specifications;
};

