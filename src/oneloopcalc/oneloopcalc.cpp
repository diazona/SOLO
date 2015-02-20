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
#include <bitset>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <openssl/sha.h>
#include "git_revision.h"
#include "../exceptions.h"
#include "../mstwpdf.h"
#include "../dss_pinlo/dss_pinlo.h"
#include "../coupling.h"
#include "../gluondist/gluondist.h"
#include "../configuration/context.h"
#include "../integration/integrationcontext.h"
#include "../hardfactors/hardfactors_position.h"
#include "../hardfactors/hardfactors_momentum.h"
#include "../integration/integrator.h"
#include "../utils/utils.h"
#include "../log.h"
#include "../hardfactors/hardfactors_radial.h"
#include "../hardfactors/hardfactor_group_parser.h"
#include "../hardfactors/hardfactor_parsed.h"

using namespace std;

/* The following functions are callback functions to be used with
 * Integrator.set_callback(). These would be invoked every time
 * the Monte Carlo routine evaluates the function.
 */

namespace trace_variable {
enum {
#define process(v) v,
#include "../integration/ictx_var_list.inc"
#undef process
    COUNT
};
}

static bitset<trace_variable::COUNT> trace_vars;

/**
 * A callback function that prints out a bunch of kinematic variables.
 */
void write_data_point(const IntegrationContext* ictx, const double real, const double imag) {
    static ofstream trace_stream("trace.output");
    if (ictx == NULL) {
        trace_stream << endl;
        return;
    }
#define process(v) if (trace_vars[(size_t)trace_variable::v]) { trace_stream << ictx->v << "\t"; }
#include "../integration/ictx_var_list.inc"
#undef process
    trace_stream << real << "\t" << imag << "\t";
    trace_stream << endl;
}

/** An IntegrationContext to store the minimum values of variables */
static IntegrationContext min_ictx(NULL, NULL);
/** An IntegrationContext to store the maximum values of variables */
static IntegrationContext max_ictx(NULL, NULL);

/** Stores a property into min_ictx and/or max_ictx if it is a min or max, respectively */
#define process(property) \
  min_ictx.property = min_ictx.property == 0 ? ictx->property : min(min_ictx.property, ictx->property); \
  max_ictx.property = max_ictx.property == 0 ? ictx->property : max(max_ictx.property, ictx->property);

/**
 * A callback function that iterates through various variables and
 * stores each into min_ictx if it is the lowest such value seen,
 * or into max_ictx if it is the highest such value seen. This is
 * used with the --minmax command-line option that allows printing
 * out the range each variable takes on during the integration.
 */
void store_minmax(const IntegrationContext* ictx, const double real, const double imag) {
    if (ictx == NULL) {
        return;
    }
#include "../integration/ictx_var_list.inc"
}
#undef process

/**
 * A callback function that writes out the result of the integration
 * if either the real or imaginary part is nonzero.
 */
void write_nonzero(const IntegrationContext* ictx, const double real, const double imag) {
    if (real != 0 || imag != 0) {
        cerr << real << "\t" << imag << endl;
    }
}

/* The following functions are callback functions to be used with
 * Integrator.set_*_callback(). These would be invoked when the
 * integration routine returns a final value, or an intermediate
 * value in the case of VEGAS, not every time it evaluates the function.
 */

/**
 * A callback for cubature integration that prints out the result of the
 * integration with its error bound.
 */
void cubature_eprint_callback(double* p_result, double* p_abserr) {
    cerr << "cubature output: " << *p_result << " err: " << *p_abserr << endl;
}
/**
 * A callback for VEGAS integration that prints out the result of the
 * integration with its error bound and chi-squared value.
 */
void vegas_eprint_callback(double* p_result, double* p_abserr, gsl_monte_vegas_state* s) {
    cerr << "VEGAS output: " << *p_result << " err: " << *p_abserr << " chisq:" << gsl_monte_vegas_chisq(s) << endl;
}
/**
 * A callback for MISER integration that prints out the result of the
 * integration with its error bound.
 */
void miser_eprint_callback(double* p_result, double* p_abserr, gsl_monte_miser_state* s) {
    cerr << "MISER output: " << *p_result << " err: " << *p_abserr << endl;
}
/**
 * A callback for quasi Monte Carlo integration that prints out the result of the
 * integration with its error bound.
 */
void quasi_eprint_callback(double* p_result, double* p_abserr, quasi_monte_state* s) {
    cerr << "QUASI output: " << *p_result << " err: " << *p_abserr << endl;
}

class ProgramConfiguration;

/**
 * Stores the results of the integration and contains methods to run the calculation.
 */
class ResultsCalculator {
private:
    /** Collection of the contexts to be used for the calculation */
    ContextCollection& cc;
    /** The thread-local context to be used for the calculation */
    ThreadLocalContext& tlctx;
    /**
     * The list of groups of hard factors
     */
    vector<const HardFactorList*> hfgroups;
    /**
     * The list of names of the hard factor groups. Each name goes with the
     * hard factor group at the corresponding index in hfgroups.
     */
    vector<string> hfgnames;
    /**
     * The list of names of the hard factors. They are stored according to the
     * group they were given in, then by order within the group.
     */
    vector<string> hfnames;

    /** The number of hard factor groups */
    size_t _hfglen;
    /** The number of hard factors */
    size_t _hflen;
    /** The length of the results arrays */
    size_t result_array_len;
    /** Flags the indices of results which have been successfully computed so far */
    bool* _valid;
    /** Array to hold the real parts of the results */
    double* real;
    /** Array to hold the imaginary parts of the results */
    double* imag;
    /** Array to hold the error bounds of the results */
    double* error;

    friend ostream& operator<<(ostream&, ResultsCalculator&);
public:
    /** Whether to trace execution */
    const bool trace;
    /** Whether to store minimum and maximum values */
    const bool minmax;
    /** Whether to calculate individual hard factors separately */
    const bool separate;

    ResultsCalculator(ContextCollection& cc, ThreadLocalContext& tlctx, ProgramConfiguration& pc);
    ~ResultsCalculator();
    /**
     * Turns a context index and a hard factor group index into an index into a
     * 1D row-major array
     */
    size_t index_from(size_t ccindex, size_t hfindex);
    /**
     * Return whether the given combination of context index and hard factor
     * group index is valid - that is, whether a result has been computed
     * for that combination
     */
    bool valid(size_t ccindex, size_t hfindex);
    /**
     * Places the result at the given context index and hard factor group
     * index into the variables real, imag, and error. This should only
     * be called after calculate().
     */
    void result(size_t ccindex, size_t hfindex, double* real, double* imag, double* error);
    /**
     * Runs the calculation.
     */
    void calculate();
private:
    /**
     * Construct an Integrator and use it
     */
    void integrate_hard_factor(const Context& ctx, const ThreadLocalContext& tlctx, const HardFactorList& hflist, size_t index);

    double xg_min, xg_max;
};

/**
 * Stores high-level program configuration variables, e.g. information
 * about which command line options were passed.
 */
class ProgramConfiguration {
public:
    ProgramConfiguration(int argc, char** argv);
    ~ProgramConfiguration();
    /**
     * Return a ContextCollection constructed using the information in
     * the ProgramConfiguration.
     */
    ContextCollection& context_collection() {
        return cc;
    }
    /**
     * Parse the hard factor specifications collected in the constructor.
     *
     * This will force a call to ContextCollection::create_contexts if
     * that has not already been done manually.
     */
    void parse_hf_specs();

    friend class ResultsCalculator;
private:
    /** Indicates whether the --trace option was specified */
    bool trace;
    /** Indicates whether the --trace-gdist option was specified */
    bool trace_gdist;
    /** Indicates whether the --minmax option was specified */
    bool minmax;
    /** Indicates whether the --separate option was specified */
    bool separate;
    /**
     * The collection of contexts to be used in the calculation. Information
     * collected from the command line options and read from configuration files
     * specified on the command line will be stored in this.
     */
    ContextCollection cc;
    /**
     * The list of transverse momenta given on the command line, if any
     */
    vector<string> pT;
    /**
     * The list of hard factor groups given on the command line
     */
    vector<const HardFactorList*> hfgroups;
    /**
     * The list of names of hard factor groups given on the command line
     */
    vector<string> hfgnames;
    /**
     * The list of names of hard factors from the groups given on the command line
     */
    vector<string> hfnames;
    vector<string> hfspecs;

    double xg_min, xg_max;
};

ProgramConfiguration::ProgramConfiguration(int argc, char** argv) : trace(false), trace_gdist(false), minmax(false), separate(false), xg_min(0), xg_max(1) {
    string gdist_type;
    bool current_arg_is_config_line = false;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (current_arg_is_config_line) {
            cc.read_config_line(a);
            current_arg_is_config_line = false;
        }
        else if (a.compare(0, 10, "--ygrange=") == 0) {
            vector<string> v = split(a, "=", 2);
            if (v.size() == 2) {
                vector<string> r = split(v[1], ":", 2);
                if (r.size() == 2) {
                    xg_min = exp(-atof(r[1].c_str()));
                    xg_max = exp(-atof(r[0].c_str()));
                    if (xg_min > xg_max) {
                        cerr << "WARNING: reversing inverted range for ln(1/xg)" << endl;
                        double t = xg_min;
                        xg_min = xg_max;
                        xg_max = t;
                    }
                }
                else {
                    cerr << "invalid range for ln(1/xg): " << v[1] << endl;
                }
            }
        }
        else if (a.compare(0, 8, "--trace=") == 0) {
            vector<string> v = split(a, "=", 2);
            if (v.size() == 2) {
                if (v[1] == "*" || v[1] == "all") {
                    trace_vars.set();
                }
                else {
                    vector<string> w = split(v[1], ",");
                    for (vector<string>::iterator it = w.begin(); it != w.end(); it++) {
                        bool handled = false;
#define process(v) if (*it == #v) { trace_vars.set((size_t)trace_variable::v); handled = true; }
#include "../integration/ictx_var_list.inc"
#undef process
                        if (!handled) {
                            cerr << "unknown trace variable " << *it << endl;
                        }
                    }
                }
                trace = trace_vars.any();
            }
        }
        else if (a == "--trace") {
            trace = true;
        }
        else if (a == "--trace-gdist") {
            trace_gdist = true;
        }
        else if (a == "--minmax") {
            minmax = true;
        }
        else if (a == "--separate") {
            separate = true;
        }
        else if (a == "-o" || a == "--option") {
            current_arg_is_config_line = true;
        }
        else if (a.compare(0, 2, "-o") == 0) {
            string s(a.substr(2));
            cc.read_config_line(s);
        }
        else if (a.compare(0, 8, "--option") == 0) {
            string s(a.substr(8));
            cc.read_config_line(s);
        }
        else if (a == "MV" || a == "fMV" ||  a == "GBW") {
            gdist_type = a;
        }
        else if (::isdigit(a[0])) {
            vector<string> pTnums = split(a, ",");
            for (vector<string>::iterator it = pTnums.begin(); it != pTnums.end(); it++) {
                pT.push_back(trim(*it, " \t"));
            }
        }
        else {
            // try opening as a file
            ifstream config;
            config.open(a.c_str());
            if (config.good()) {
                logger << "Reading config file " << a << endl;
                config >> cc;
                config.close();
            }
            else {
                hfspecs.push_back(a);
            }
        }
    }

    if (!pT.empty()) {
        cc.erase("pT");
        for (vector<string>::iterator it = pT.begin(); it != pT.end(); it++) {
            cc.add("pT", *it);
        }
    }
    if (!gdist_type.empty()) {
        cc.set("gdist", gdist_type);
    }
    cc.trace_gdist = trace_gdist;
    if (hfspecs.empty()) {
        hfspecs.push_back("lo");
        hfspecs.push_back("nlo");
    }
}

void ProgramConfiguration::parse_hf_specs() {
    if (cc[0].hardfactor_definitions.empty()) {
        throw MissingPropertyException("no hard factors defined");
    }

    HardFactorParser parser;
    for (vector<string>::const_iterator it = cc[0].hardfactor_definitions.begin(); it != cc[0].hardfactor_definitions.end(); it++) {
        parser.parse_file(*it);
    }

    assert(hfgroups.empty());
    for (vector<string>::const_iterator it = hfspecs.begin(); it!= hfspecs.end(); it++) {
        const ParsedHardFactorGroup* hfg = ParsedHardFactorGroup::parse(*it, cc[0].exact_kinematics);
        hfgroups.push_back(hfg->objects);
        hfgnames.push_back(hfg->label);
        hfnames.insert(hfnames.end(), hfg->specifications.begin(), hfg->specifications.end());
    }
    assert(!hfgroups.empty());
    assert(hfgroups.size() == hfgnames.size());
    assert(hfnames.size() >= hfgnames.size());
}

ProgramConfiguration::~ProgramConfiguration() {
    for (vector<const HardFactorList*>::iterator it = hfgroups.begin(); it != hfgroups.end(); it++) {
        delete *it;
    }
}

ResultsCalculator::ResultsCalculator(ContextCollection& cc, ThreadLocalContext& tlctx, ProgramConfiguration& pc) :
    cc(cc), tlctx(tlctx), hfgroups(pc.hfgroups), hfgnames(pc.hfgnames), hfnames(pc.hfnames), result_array_len(cc.size()), _hfglen(0), _hflen(0), trace(pc.trace), minmax(pc.minmax), separate(pc.separate), xg_min(pc.xg_min), xg_max(pc.xg_max) {
    _hfglen = hfgroups.size();
    assert(_hfglen > 0);
    if (separate) {
        for (vector<const HardFactorList*>::iterator hflit = hfgroups.begin(); hflit != hfgroups.end(); hflit++) {
            _hflen += (*hflit)->size();
        }
        result_array_len *= _hflen;
    }
    else {
        result_array_len *= _hfglen;
    }
    _valid = new bool[result_array_len];
    real = new double[result_array_len];
    imag = new double[result_array_len];
    error = new double[result_array_len];
    fill(_valid, _valid + result_array_len, false);
}
ResultsCalculator::~ResultsCalculator() {
    delete[] _valid;
    delete[] real;
    delete[] imag;
    delete[] error;
}
size_t ResultsCalculator::index_from(size_t ccindex, size_t hfindex) {
    size_t index = ccindex * (separate ? _hflen : _hfglen) + hfindex;
    assert(index < result_array_len);
    return index;
}
bool ResultsCalculator::valid(size_t ccindex, size_t hfindex) {
    return _valid[index_from(ccindex, hfindex)];
}

void ResultsCalculator::result(size_t ccindex, size_t hfindex, double* real, double* imag, double* error) {
    size_t index = index_from(ccindex, hfindex);
    if (valid(ccindex, hfindex)) {
        *real = this->real[index];
        *imag = this->imag[index];
        *error = this->error[index];
    }
    else {
        ostringstream s;
        s << "Invalid results at ccindex " << ccindex << ", hfindex " << hfindex << endl;
        throw s.str().c_str();
    }
}

void ResultsCalculator::calculate() {
    size_t cc_index = 0, hf_index = 0;
    for (ContextCollection::iterator it = cc.begin(); it != cc.end(); it++) {
        Context ctx = *it;
        cerr << "Beginning calculation at pT = " << sqrt(ctx.pT2) << ", Y = " << ctx.Y << endl;
        try {
            hf_index = 0;
            // recall the definition
            // typedef HardFactorList std::vector<const HardFactor*>
            for (vector<const HardFactorList*>::iterator hit = hfgroups.begin(); hit != hfgroups.end(); hit++) {
                if (separate) {
                    // go through the hard factors in each group one at a time
                    HardFactorList one_hf;
                    for (HardFactorList::const_iterator hfit = (*hit)->begin(); hfit != (*hit)->end(); hfit++) {
                        one_hf.assign(1, *hfit);
                        integrate_hard_factor(ctx, tlctx, one_hf, index_from(cc_index, hf_index));
                        hf_index++;
                    }
                }
                else {
                    integrate_hard_factor(ctx, tlctx, **hit, index_from(cc_index, hf_index));
                    hf_index++;
                }
            }
            cerr << "...done" << endl;
        }
        catch (const exception& e) {
            cerr << e.what() << endl;
        }
        cc_index++;
    }
}

void ResultsCalculator::integrate_hard_factor(const Context& ctx, const ThreadLocalContext& tlctx, const HardFactorList& hflist, size_t index) {
    double l_real, l_imag, l_error;
    Integrator integrator(&ctx, &tlctx, hflist, xg_min, xg_max);
    if (trace) {
        integrator.set_callback(write_data_point);
    }
    else if (minmax) {
        integrator.set_callback(store_minmax);
    }
    integrator.set_cubature_callback(cubature_eprint_callback);
    integrator.set_miser_callback(miser_eprint_callback);
    integrator.set_vegas_callback(vegas_eprint_callback);
    integrator.set_quasi_callback(quasi_eprint_callback);
    integrator.integrate(&l_real, &l_imag, &l_error);
    real[index] = l_real;
    imag[index] = l_imag;
    error[index] = l_error;
    _valid[index] = true;
}

/**
 * Write the list of results in a ResultsCalculator to the given output stream.
 */
ostream& operator<<(ostream& out, ResultsCalculator& rc) {
    using std::setw;
    const string OFS = " ";   // output field separator - make sure fields are space-separated
    const string BLANK = " "; // a blank field
    size_t lw = 6;
    size_t rw = 14; // "label width" and "result width"

    // write headers
    if (rc.separate) {
        out << setw(lw) << left << "pT" << OFS << setw(lw) << "Y" << OFS;
        for (size_t hfgindex = 0; hfgindex < rc._hfglen; hfgindex++) {
            out << setw(rw) << rc.hfgnames[hfgindex] << OFS;
            size_t hflen = rc.hfgroups[hfgindex]->size();
            for (size_t hfindex = 1; hfindex < 2 * hflen; hfindex++) {
                out << setw(rw) << BLANK << OFS;
            }
        }
        out << setw(rw) << "total" << endl;

        out << setw(lw) << BLANK << OFS << setw(lw) << BLANK << OFS;
        for (vector<string>::iterator termname_iterator = rc.hfnames.begin(); termname_iterator != rc.hfnames.end(); termname_iterator++) {
            ostringstream valstream;
            valstream << *termname_iterator << "-val";
            out << setw(rw) << valstream.str() << OFS;
            ostringstream errstream;
            errstream << *termname_iterator << "-err";
            out << setw(rw) << errstream.str() << OFS;
        }
        out << endl;
    }
    else {
        out << setw(lw) << left << "pT" << OFS << setw(lw) << "Y" << OFS;
        for (vector<string>::iterator it = rc.hfgnames.begin(); it != rc.hfgnames.end(); it++) {
            ostringstream valstream;
            valstream << *it << "-val";
            out << setw(rw) << valstream.str() << OFS;
            ostringstream errstream;
            errstream << *it << "-err";
            out << setw(rw) << errstream.str() << OFS;
        }
        out << setw(rw) << "total" << endl;
    }

    // write data
    double l_real, l_imag, l_error;
    bool all_valid = true;
    for (size_t ccindex = 0; ccindex < rc.cc.size(); ccindex++) {
        out << setw(lw) << sqrt(rc.cc[ccindex].pT2) << OFS;
        out << setw(lw) << rc.cc[ccindex].Y << OFS;

        double total = 0;
        size_t hfglen = rc.separate ? rc._hflen : rc._hfglen;
        bool row_valid = true;
        for (size_t hfgindex = 0; hfgindex < hfglen; hfgindex++) {
            if (rc.valid(ccindex, hfgindex)) {
                rc.result(ccindex, hfgindex, &l_real, &l_imag, &l_error);
                out << setw(rw) << l_real << OFS << setw(rw) << l_error << OFS;
                total += l_real;
            }
            else {
                out << setw(rw) << "---" << OFS << setw(rw) << "---" << OFS;
                all_valid = row_valid = false;
            }
        }
        if (row_valid) {
            out << setw(rw) << total << endl;
        }
        else {
            out << setw(rw) << "---" << endl;
        }
    }
    if (!all_valid) {
        out << "WARNING: some results were not computed" << endl;
    }

    if (rc.minmax) {
#define process(v) out << #v << "\t" << min_ictx.v << "\t" << max_ictx.v << "\t" << endl;
#include "../integration/ictx_var_list.inc"
#undef process
    }
    return out;
}

/* The output stream to write logging messages to. Declared in log.h. */
ostream& logger = cerr;

/** The one instance of ResultsCalculator used for the program */
static ResultsCalculator* p_rc = NULL;

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
        os << setw(2) << (unsigned int)*ptr;
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
        SHA1_Update(&c, buffer, i.gcount());
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
    ContextCollection cc = pc.context_collection();
    if (cc.empty()) {
        logger << "No momenta or no rapidities specified!" << endl;
        return 1;
    }

    /* First write out all the configuration variables. Having the configuration written
     * out as part of the output file makes it easy to tell what parameters were used in
     * and given run, and is also useful in case we want to reproduce a run.
     */
#ifdef GIT_REVISION
    cout << "# git revision " << GIT_REVISION;
#ifdef GIT_DIRTY
    cout << " (dirty)";
#endif
    cout << endl;
#endif
    {
        FileDataGluonDistribution* fgdist = dynamic_cast<FileDataGluonDistribution*>(cc[0].gdist);
        if (fgdist != NULL) {
            // print hashes of the input files
            // TODO: make the gdist compute the hashes itself
            cout << "# momentum gdist file hash: " << sha1_file(cc.get("gdist_momentum_filename")) << endl;
            cout << "# position gdist file hash: " << sha1_file(cc.get("gdist_position_filename")) << endl;
        }
        string hf_definition_filename = cc.get("hf_definitions");
        if (!hf_definition_filename.empty()) {
            cout << "# hard factor definition file hash: " << sha1_file(hf_definition_filename) << endl;
            ifstream hfdefs(hf_definition_filename.c_str());
            if (!hfdefs) {
                ostringstream oss;
                oss << "Error opening hard factor definition file: " << hf_definition_filename;
                throw ios_base::failure(oss.str());
            }
            cerr << "BEGIN hf definition file" << endl << hfdefs.rdbuf() << "END hf definition file" << endl;
        }
    }

#ifdef EXACT_LIMIT_SCHEME
    cout << "# EXACT_LIMIT_SCHEME = " << EXACT_LIMIT_SCHEME << endl;
#else
    cout << "# EXACT_LIMIT_SCHEME undefined" << endl;
#endif

    cout << cc << "------------" << endl;

    cc.create_contexts();
    if (cc.empty()) {
        logger << "No valid momentum/rapidity combinations specified!" << endl;
        return 1;
    }

    // parse hard factor specifications after creating contexts
    pc.parse_hf_specs();

    // Only create ThreadLocalContext here because cc may not have
    // values for pdf_filename and ff_filename before create_contexts
    // is called
    ThreadLocalContext tlctx(cc);

    // pc.parse_hf_specs() needs to happen before this
    ResultsCalculator rc(cc, tlctx, pc);
    p_rc = &rc;

    /* Set up a signal handler so that if the program receives a SIGINT (Ctrl+C)
     * or SIGTERM (e.g. runs out of time in PBS), it will invoke termination_handler()
     * to print what results it has so far
     */
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
        string spaces(e.GetPos(), ' ');
        cerr << "Parser error: " << e.GetMsg() << endl;
        cerr << "in expression:" << endl;
        cerr << e.GetExpr() << endl;
        cerr << spaces << "^" << endl;
        return 1;
    }
    catch (const exception& e) {
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
