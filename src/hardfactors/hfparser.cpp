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

#include <iostream>
#include "hardfactor.h"
#include "hardfactor_parsed.h"
#include "hardfactor_group_parser.h"

using std::vector;

HardFactorTermList parse_hf_definitions(const char* filename) {
    HardFactorParser parser;
    HardFactorTermList h;
    try {
        h = parser.parse(filename);
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << endl;
        // h will be empty at this point
        assert(h.empty());
    }
    return h;
}

int main(const int argc, char** argv) {
    if (argc < 2) {
        return 2;
    }
    HardFactorTermList hl = parse_hf_definitions(argv[1]);
    if (hl.empty()) {
        return 0;
    }
    else {
        std::cerr << "Error parsing definitions" << endl;
        return 1;
    }
}