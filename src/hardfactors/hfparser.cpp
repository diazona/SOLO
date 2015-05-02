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
#include <typeinfo>
#include <vector>
#include <muParser.h>
#include "hardfactor.h"
#include "hardfactor_parser.h"

using std::vector;
using std::cout;
using std::cerr;
using std::endl;

static bool encountered_error;
static bool debug;

static HardFactorList hl;
static std::vector<const HardFactorGroup*> hgl;

bool print_err(const std::exception& e, const std::string& filename, const size_t line_number) {
    cerr << "Error parsing at line " <<  line_number << " of " << filename << endl;
    const InvalidHardFactorDefinitionException* ihfde = dynamic_cast<const InvalidHardFactorDefinitionException*>(&e);
    if (ihfde != NULL) {
        cerr << ihfde->line << endl;
    }
    cerr << e.what() << endl;
    encountered_error = true;
    return false;
}

void handle_hard_factor(const HardFactor& hf) {
    if (debug) {
        cerr << "parsed hard factor " << hf.get_name() << endl;
    }
    hl.push_back(&hf);
}

void handle_hard_factor_group(const HardFactorGroup& hfg) {
    if (debug) {
        cerr << "parsed hard factor group " << hfg.label << endl;
    }
    hgl.push_back(&hfg);
}

int main(const int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " filename [filename...]" << endl;
        return 2;
    }
    HardFactorRegistry registry;
    HardFactorParser parser(registry);
    vector<string> files;
    bool verbose = false;
    encountered_error = false;
    debug = false;
    for (size_t i = 1; i < static_cast<size_t>(argc); i++) {
        string s(argv[i]);
        if (s.compare(0, 2, "--") == 0) {
            if (s == "--debug") {
            }
            else if (s == "--verbose") {
                verbose = true;
            }
        }
        else {
            files.push_back(s);
        }
    }
    parser.set_error_handler(print_err);
    parser.set_hard_factor_callback(handle_hard_factor);
    parser.set_hard_factor_group_callback(handle_hard_factor_group);
    for (vector<string>::const_iterator it = files.begin(); it != files.end(); it++) {
        try {
            parser.parse_file(*it);
        }
        catch (const std::exception& e) {
            cerr << e.what() << endl;
            return 1;
        }
        catch (const mu::ParserError& e) {
            cerr << e.GetMsg() << endl;
            return 1;
        }
    }
    parser.flush_groups();
    if (encountered_error) {
        return 1;
    }
    cout << "Found " << hl.size() << " hard factors and " << hgl.size() << " groups:" << endl;
    size_t wrap = 0;
    for (HardFactorList::const_iterator it = hl.begin(); it != hl.end(); it++) {
        const HardFactor& hf = **it;
        if (verbose) {
            cout << "- " << hf.get_name() << endl << "  ";
            switch (hf.get_order()) {
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
            try {
                const ParsedHardFactorTerm& phft = dynamic_cast<const ParsedHardFactorTerm&>(hf);
                cout << "  Fs = " << phft.Fs_expr() << endl;
                cout << "  Fn = " << phft.Fn_expr() << endl;
                cout << "  Fd = " << phft.Fd_expr() << endl;
                cout << endl;
            }
            catch (const bad_cast& e) {
                try {
                    const ParsedCompositeHardFactor& pchf = dynamic_cast<const ParsedCompositeHardFactor&>(hf);
                    cout << "  contains ";
                    const size_t nterms = pchf.get_term_count();
                    const HardFactorTerm* const* hfl = pchf.get_terms();
                    for (size_t i = 0; i < nterms; i++) {
                        if (i > 0) {
                            cout << ", ";
                        }
                        cout << hfl[i]->get_name();
                    }
                    cout << endl << endl;
                }
                catch (const bad_cast& e) {
                    cerr << "Unknown hard factor type" << endl;
                }
            }
        }
        else {
            if (it != hl.begin()) {
                cout << ", ";
                if (wrap % 6 == 0) {
                    cout << endl;
                }
            }
            cout << hf.get_name();
            wrap++;
        }
    }
    if (!verbose) {
        cout << endl;
    }
    for (std::vector<const HardFactorGroup*>::const_iterator it = hgl.begin(); it != hgl.end(); it++) {
        const HardFactorGroup& hfg = **it;
        cout << "+ " << hfg.label << " : ";
        for (HardFactorList::const_iterator lit = hfg.objects->begin(); lit != hfg.objects->end(); lit++) {
            const HardFactor& hf = **lit;
            if (lit != hfg.objects->begin()) {
                cout << ", ";
            }
            cout << hf.get_name();
        }
        cout << endl;
    }
    return 0;
}