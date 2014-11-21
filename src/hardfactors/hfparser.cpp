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

bool print_err(const std::exception& e, const std::string& filename, const size_t line_number) {
    std::cerr << "Error parsing at line " <<  line_number << " of " << filename << endl;
    const InvalidHardFactorDefinitionException* ihfde = dynamic_cast<const InvalidHardFactorDefinitionException*>(&e);
    if (ihfde != NULL) {
        std::cerr << ihfde->line << endl;
    }
    std::cerr << e.what() << endl;
    return false;
}

HardFactorTermList parse_hf_definitions(const char* filename) {
    HardFactorParser parser;
    try {
        parser.parse_file(filename);
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << endl;
    }
    return parser.get_hard_factors();
}

int main(const int argc, char** argv) {
    if (argc < 2) {
        return 2;
    }
    HardFactorTermList hl = parse_hf_definitions(argv[1]);
    if (hl.empty()) {
        std::cerr << "No definitions found";
        return 0;
    }
    else {
        return 0;
    }
}