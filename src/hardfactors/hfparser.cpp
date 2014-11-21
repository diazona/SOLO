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
#include <muParser.h>
#include "hardfactor.h"
#include "hardfactor_parsed.h"
#include "hardfactor_group_parser.h"

using std::vector;
using std::cerr;
using std::endl;

bool print_err(const std::exception& e, const std::string& filename, const size_t line_number) {
    cerr << "Error parsing at line " <<  line_number << " of " << filename << endl;
    const InvalidHardFactorDefinitionException* ihfde = dynamic_cast<const InvalidHardFactorDefinitionException*>(&e);
    if (ihfde != NULL) {
        cerr << ihfde->line << endl;
    }
    cerr << e.what() << endl;
    return false;
}

int main(const int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " filename [filename...]" << endl;
        return 2;
    }
    HardFactorParser parser;
    for (size_t i = 1; i < argc; i++) {
        try {
            parser.parse_file(argv[i], print_err);
        }
        catch (const std::exception& e) {
            cerr << e.what() << endl;
        }
        catch (const mu::ParserError& e) {
            cerr << e.GetMsg() << endl;
        }
    }
    HardFactorTermList hl = parser.get_hard_factors();
    cout << "Found " << hl.size() << " hard factors:" << endl;
    for (HardFactorTermList::const_iterator it = hl.begin(); it != hl.end(); it++) {
        cout << "  " << (*it)->get_name() << endl;
    }
    return 0;
}