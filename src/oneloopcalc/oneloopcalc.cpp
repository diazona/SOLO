/*
 * A calculation of the NLO cross section of pA->pion collisions
 *
 * This file contains the driver code, including main() and some other stuff
 *
 * Copyright 2012 David Zaslavsky
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

#include <cassert>
#include <cctype>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <openssl/sha.h>
#include "muParserError.h"
#include "git_revision.h"
#include "../exceptions.h"
#include "../configuration/context.h"
#include "../utils/utils.h"
#include "../log.h"
#include "programconfiguration.h"
#include "resultscalculator.h"

using namespace std;

/* The output stream to write logging messages to. Declared in log.h. */
ostream& logger = cerr;

/**
 * The one instance of ResultsCalculator used for the program.
 *
 * This has to be a pointer, not a reference, because it's undefined
 * until the actual ResultsCalculator object is constructed in the
 * run() function.
 */
static ResultsCalculator* p_rc;

/**
 * Takes care of finishing the program if it gets interrupted by a signal.
 *
 * This happens when a PBS job is cut off before it finishes, for example. This
 * function will write out all results computed so far by writing the
 * ResultsCalculator object to standard output, and then exit the program.
 */
void termination_handler(int signal) {
    static bool terminated = false;
    if (!terminated) {
        terminated = true;
        cout << *p_rc;

        time_t rawtime;
        time(&rawtime);
        logger << "Terminating at " << ctime(&rawtime) << endl;
    }
    exit(2);
}

// from http://stackoverflow.com/questions/3969047/is-there-a-standard-way-of-representing-an-sha1-hash-as-a-c-string-and-how-do-i
string get_hex_representation(const unsigned char* bytes, size_t length) {
    ostringstream os;
    os.fill('0');
    os << hex;
    for(const unsigned char * ptr=bytes; ptr < bytes+length; ptr++) {
        os << setw(2) << static_cast<unsigned int>(*ptr);
    }
    return os.str();
}

string sha1_file(string filename) {
    char buffer[1024];
    unsigned char hash[SHA_DIGEST_LENGTH];
    SHA_CTX c;
    SHA1_Init(&c);
    ifstream i(filename.c_str());
    if (!i) {
        ostringstream oss;
        oss << "Error opening file for SHA checksum: " << filename;
        throw ios_base::failure(oss.str());
    }
    while (i) {
        i.read(buffer, sizeof(buffer));
        SHA1_Update(&c, buffer, static_cast<size_t>(i.gcount()));
    }
    i.close();
    SHA1_Final(hash, &c);
    return get_hex_representation(hash, SHA_DIGEST_LENGTH);
}

/**
 * GSL error handler function that throws a GSLException.
 */
void gsl_error_throw(const char* reason, const char* file, int line, int gsl_errno) {
    throw GSLException(reason, file, line, gsl_errno);
}

/**
 * Runs the program.
 *
 * This is like main() except that it can throw exceptions, which will
 * be caught in the real main().
 */
int run(int argc, char** argv) {
    time_t rawtime;

    time(&rawtime);
    logger << "Starting at " << ctime(&rawtime) << endl;

    gsl_set_error_handler(&gsl_error_throw);

    ProgramConfiguration pc(argc, argv);

    /* First write out all the configuration variables. Having the configuration written
     * out as part of the output file makes it easy to tell what parameters were used in
     * and given run, and is also useful in case we want to reproduce a run.
     */
    if (pc.print_config()) {
#ifdef GIT_REVISION
        cout << "# git revision " << GIT_REVISION;
#ifdef GIT_DIRTY
        cout << " (dirty)";
#endif
        cout << endl;
#endif
    }

    ResultsCalculator rc(pc);

    {
        FileDataGluonDistribution* fgdist = dynamic_cast<FileDataGluonDistribution*>(rc.cc[0].gdist);
        if (pc.print_config() && fgdist != NULL) {
            // print hashes of the input files
            // TODO: make the gdist compute the hashes itself
            cout << "# momentum gdist file hash: " << sha1_file(rc.cc.config().get("gdist_momentum_filename")) << endl;
            cout << "# position gdist file hash: " << sha1_file(rc.cc.config().get("gdist_position_filename")) << endl;
        }

        vector<string> hfdefs = rc.cc[0].hardfactor_definitions;
        for (vector<string>::const_iterator it = hfdefs.begin(); it != hfdefs.end(); it++) {
            string hf_definition_filename = *it;
            ifstream hfdefs(hf_definition_filename.c_str());
            if (!hfdefs) {
                ostringstream oss;
                oss << "Error opening hard factor definition file: " << hf_definition_filename;
                throw ios_base::failure(oss.str());
            }
            cerr << "BEGIN hf definition file " << hf_definition_filename << endl << hfdefs.rdbuf() << "END hf definition file " << hf_definition_filename << endl;
            hfdefs.close();
            if (pc.print_config()) {
                cout << "# hard factor definition file hash: " << hf_definition_filename << ": " << sha1_file(hf_definition_filename) << endl;
            }
        }
    }

    if (pc.print_config()) {
#ifdef EXACT_LIMIT_SCHEME
        cout << "# EXACT_LIMIT_SCHEME = " << EXACT_LIMIT_SCHEME << endl;
#else
        cout << "# EXACT_LIMIT_SCHEME undefined" << endl;
#endif
    }

    if (pc.print_config()) {
        cout << rc.cc.config() << "------------" << endl;
    }

    if (rc.cc.empty()) {
        logger << "No valid momentum/rapidity combinations specified!" << endl;
        return 1;
    }

    /* Set up a signal handler so that if the program receives a SIGINT (Ctrl+C)
     * or SIGTERM (e.g. runs out of time in PBS), it will invoke termination_handler()
     * to print what results it has so far
     */
    p_rc = &rc;
    struct sigaction siga;
    siga.sa_handler = termination_handler;
    struct sigaction oldsiga;
    sigaction(SIGTERM, &siga, &oldsiga);
    sigaction(SIGINT, &siga, &oldsiga);
    // Run the actual calculation
    rc.calculate();
    // Reset the signal handler
    sigaction(SIGTERM, &oldsiga, NULL);
    sigaction(SIGINT, &oldsiga, NULL);
    // And print out results
    cout << rc;

    time(&rawtime);
    logger << "Ending at " << ctime(&rawtime) << endl;

    return 0;
}

/**
 * This just calls run() and catches any exceptions that may be thrown.
 */
int main(int argc, char** argv) {
    try {
        return run(argc, argv);
    }
    catch (const mu::ParserError& e) {
        cerr << "Parser error: " << e.GetMsg() << endl;
        string expr = e.GetExpr();
        if (!expr.empty()) {
            cerr << "in expression:" << endl;
            cerr << expr << endl;
            string spaces(e.GetPos(), ' ');
            cerr << spaces << "^" << endl;
        }
        return 1;
    }
    catch (const std::exception& e) {
        cerr << "Caught exception:" << endl << e.what() << endl;
        return 1;
    }
    catch (const char* c) {
        cerr << "Caught error message:" << endl << c << endl;
        return 1;
    }
    catch (...) {
        cerr << "Unknown error" << endl;
        return 1;
    }
}
