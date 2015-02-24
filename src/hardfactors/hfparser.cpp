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
    HardFactorList hl = parser.get_hard_factors();
    cout << "Found " << hl.size() << " hard factors:" << endl;
    for (HardFactorList::const_iterator it = hl.begin(); it != hl.end(); it++) {
        const HardFactor* hf = *it;
        if (hf != NULL) {
            cout << "- " << hf->get_name() << endl << "  ";
            switch ((*it)->get_order()) {
                case HardFactor::LO:
                    cout << "LO";
                    break;
                case HardFactor::NLO:
                    cout << "NLO";
                    break;
                case HardFactor::MIXED:
                    cout << "mixed";
                    break;
                default:
                    cout << "unknown!";
                    break;
            }
            cout << endl;
        }
        const ParsedHardFactorTerm* phft = dynamic_cast<const ParsedHardFactorTerm*>(hf);
        if (phft != NULL) {
            cout << "  Fs = " << phft->Fs_expr() << endl;
            cout << "  Fn = " << phft->Fn_expr() << endl;
            cout << "  Fd = " << phft->Fd_expr() << endl;
            cout << endl;
        }
        else {
            const ParsedCompositeHardFactor* pchf = dynamic_cast<const ParsedCompositeHardFactor*>(hf);
            if (pchf != NULL) {
                cout << "  contains ";
                const size_t nterms = pchf->get_term_count();
                const HardFactorTerm* const* hfl = pchf->get_terms();
                for (size_t i = 0; i < nterms; i++) {
                    if (i > 0) {
                        cout << ", ";
                    }
                    cout << hfl[i]->get_name();
                }
                cout << endl << endl;
            }
        }
    }
    return 0;
}