/*
 * Part of SOLO
 *
 * Copyright 2012-2015 David Zaslavsky
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <iostream>
#include "../log.h"
#include "configuration.h"
#include "context.h"

ostream& logger = cerr;

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <filename.cfg> ..." << endl;
        return 1;
    }
    Configuration conf(canonicalize);
    for (size_t i = 1; i < static_cast<size_t>(argc); i++) {
        ifstream config;
        config.open(argv[i]);
        if (config.good()) {
            config >> conf;
        }
        config.close();
    }
    try {
        ContextCollection cc(conf);
        cout << "Successfully parsed files";
        if (cc.empty()) {
            cout << endl << "No contexts defined!" << endl << cc.config() << endl;
        }
        else {
            cout << " into " << cc.size() << " contexts" << endl;
            cout << cc.config() << endl;
            for (ContextCollection::iterator it = cc.begin(); it != cc.end(); it++) {
                cout << "Context: pT2 = " << it->pT2 << ", Y = " << it->Y << ", seed = " << it->pseudorandom_generator_seed << endl;
            }
        }
    }
    catch (const exception& e) {
        cout << "Error in parsing: " << e.what() << endl;
        return 1;
    }
    return 0;
}
