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
#include <iostream>
#include <iomanip>
#include <sstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include "git_revision.h"
#include "mstwpdf.h"
#include "dss_pinlo.h"
#include "coupling.h"
#include "gluondist.h"
#include "context.h"
#include "integrationcontext.h"
#include "hardfactors_position.h"
#include "hardfactors_momentum.h"
#include "integrator.h"
#include "utils.h"
#include "log.h"
#include "hardfactors_radial.h"

using namespace std;

/* The following functions are callback functions to be used with
 * Integrator.set_callback(). These would be invoked every time
 * the Monte Carlo routine evaluates the function.
 */

/**
 * A callback function that prints out a bunch of kinematic variables.
 */
void write_data_point(const IntegrationContext* ictx, const double real, const double imag) {
    if (ictx) {
//         if ((++count) % 500 == 0) {
            cerr
            << ictx->z << "\t"
            << ictx->xi << "\t"
            << ictx->xx << "\t"
            << ictx->xy << "\t"
            << ictx->yx << "\t"
            << ictx->yy << "\t"
            << ictx->bx << "\t"
            << ictx->by << "\t"
            << ictx->kT2 << "\t"
            << ictx->Qs2 << "\t"
            << ictx->xp << "\t"
            << ictx->xg << "\t"
            << ictx->qqfactor << "\t"
            << ictx->gqfactor << "\t"
            << real << "\t"
            << imag << endl;
//         }
    }
    else {
        cerr << endl;
    }
}

/** An IntegrationContext to store the minimum values of variables */
static IntegrationContext min_ictx(NULL, NULL);
/** An IntegrationContext to store the maximum values of variables */
static IntegrationContext max_ictx(NULL, NULL);

/** Stores a property into min_ictx and/or max_ictx if it is a min or max, respectively */
#define store(property) \
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
    store(z);
    store(xi);
    store(xx);
    store(xy);
    store(yx);
    store(yy);
    store(bx);
    store(by);
    store(q1x);
    store(q1y);
    store(q2x);
    store(q2y);
    store(q3x);
    store(q3y);
    store(kT);
    store(kT2);
    store(xp);
    store(xg);
    store(xiprime);
    store(Qs2);
    store(alphas);
}

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
    vector<HardFactorList*> hfgroups;
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
    /** The number of results that have been computed so far */
    size_t _valid;
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
    void integrate_hard_factor(const Context& ctx, const ThreadLocalContext& tlctx, HardFactorList& hflist, double* l_real, double* l_imag, double* l_error);
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

    friend class ResultsCalculator;
private:
    /** Indicates whether the --trace option was specified */
    bool trace;
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
    vector<HardFactorList*> hfgroups;
    /**
     * The list of names of hard factor groups given on the command line
     */
    vector<string> hfgnames;
    /**
     * The list of names of hard factors from the groups given on the command line
     */
    vector<string> hfnames;
    
    /**
     * Parse a string specification of a hard factor group and add the
     * resulting hard factors and names to the lists
     */
    void parse_hf_spec(const string& spec);
};

ProgramConfiguration::ProgramConfiguration(int argc, char** argv) : trace(false), minmax(false), separate(false) {
    string gdist_type;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a == "--trace") {
            trace = true;
        }
        else if (a == "--minmax") {
            minmax = true;
        }
        else if (a == "--separate") {
            separate = true;
        }
        else if (a == "MV" || a == "fMV" ||  a == "GBW") {
            gdist_type = a;
        }
        else if (a[0] == 'h' || a[1] == ':' || a == "lo" || a == "nlo") {
            parse_hf_spec(a);
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
                logger << "Unrecognized argument " << a << endl;
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
    if (hfgroups.empty()) {
        parse_hf_spec("lo");
        parse_hf_spec("nlo");
    }
    assert(!hfgroups.empty());
    assert(hfgroups.size() == hfgnames.size());
    assert(hfnames.size() >= hfgnames.size());
}

ProgramConfiguration::~ProgramConfiguration() {
    for (vector<HardFactorList*>::iterator it = hfgroups.begin(); it != hfgroups.end(); it++) {
        delete *it;
    }
}

/**
 * This defines the hard factor group that is used when "lo" is given
 * on the command line
 */
static const char* default_lo_spec = "m:h02qq,m:h02gg";
/**
 * This defines the hard factor group that is used when "nlo" is given
 * on the command line
 */
static const char* default_nlo_spec = "r:h12qq,m:h14qq,r:h12gg,m:h12qqbar,m:h16gg,r:h112gq,r:h122gq,m:h14gq,r:h112qg,r:h122qg,m:h14qg";

void ProgramConfiguration::parse_hf_spec(const string& spec) {
    vector<string> splitspec;
    // Split the specification string on commas to get individual hard factor names
    if (spec == "lo") {
        splitspec = split(default_lo_spec, ", ");
    }
    else if (spec == "nlo") {
        splitspec = split(default_nlo_spec, ", ");
    }
    else {
        splitspec = split(spec, ", ");
    }
    HardFactorList* hfobjs = new HardFactorList();
    assert(hfobjs != NULL);
    // Iterate over the individual hard factor names
    for (vector<string>::iterator it = splitspec.begin(); it != splitspec.end(); it++) {
        string orig_s = *it;
        string s = orig_s;
        // The default is to create a position-space hard factor
        const HardFactorRegistry* registry = position::registry::get_instance();
        if (s[1] == ':') {
            switch (s[0]) {
                case 'm':
                    registry = momentum::registry::get_instance();
                    break;
                case 'r':
                    registry = radial::registry::get_instance();
                    break;
                case 'p':
                    registry = position::registry::get_instance();
                    break;
            }
            // chop off prefix
            s = s.substr(2);
        }
        // Pass the remaining name (e.g. "h02qq") to the hard factor registry
        // to get the actual hard factor object
        const HardFactor* hf = registry->get_hard_factor(s);
        if (hf == NULL) {
            logger << "No such hard factor " << orig_s << endl;
        }
        else {
            hfobjs->push_back(hf);
        }
    }
    if (hfobjs->size() == 0) {
        logger << "No valid hard factors in spec " << spec << endl;
        throw "Error parsing hard factors";
    }
    // This constitutes one group. Add it to the list and add the
    // specification to the lists of names.
    hfgroups.push_back(hfobjs);
    hfgnames.push_back(spec);
    hfnames.insert(hfnames.end(), splitspec.begin(), splitspec.end());
}

ResultsCalculator::ResultsCalculator(ContextCollection& cc, ThreadLocalContext& tlctx, ProgramConfiguration& pc) :
    cc(cc), tlctx(tlctx), hfgroups(pc.hfgroups), hfgnames(pc.hfgnames), hfnames(pc.hfnames), _hfglen(0), _hflen(0), _valid(0), trace(pc.trace), minmax(pc.minmax), separate(pc.separate) {
    _hfglen = hfgroups.size();
    size_t cclen = cc.size();
    real = new double[_hfglen * cclen];
    imag = new double[_hfglen * cclen];
    error = new double[_hfglen * cclen];
    if (separate) {
        for (vector<HardFactorList*>::iterator hflit = hfgroups.begin(); hflit != hfgroups.end(); hflit++) {
            _hflen += (*hflit)->size();
        }
    }
}
ResultsCalculator::~ResultsCalculator() {
    delete[] real;
    delete[] imag;
    delete[] error;
}
size_t ResultsCalculator::index_from(size_t ccindex, size_t hfindex) {
    return ccindex * (separate ? _hflen : _hfglen) + hfindex;
}
bool ResultsCalculator::valid(size_t ccindex, size_t hfindex) {
    return index_from(ccindex, hfindex) < _valid;
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
        s << "Invalid index: " << index << " out of " << _valid << " results computed" << endl;
        throw s.str().c_str();
    }
}

void ResultsCalculator::calculate() {
    double* l_real = real;
    double* l_imag = imag;
    double* l_error = error;
    for (ContextCollection::iterator it = cc.begin(); it != cc.end(); it++) {
        Context ctx = *it;
        cerr << "Beginning calculation at pT = " << sqrt(ctx.pT2) << ", Y = " << ctx.Y << endl;
        try {
            // recall the definition
            // typedef HardFactorList std::vector<const HardFactor*>
            for (vector<HardFactorList*>::iterator hit = hfgroups.begin(); hit != hfgroups.end(); hit++) {
                if (separate) {
                    // go through the hard factors in each group one at a time
                    HardFactorList one_hf;
                    one_hf.resize(1);
                    for (HardFactorList::iterator hfit = (*hit)->begin(); hfit != (*hit)->end(); hfit++) {
                        one_hf[0] = *hfit;
                        integrate_hard_factor(ctx, tlctx, one_hf, l_real++, l_imag++, l_error++);
                    }
                }
                else {
                    integrate_hard_factor(ctx, tlctx, **hit, l_real++, l_imag++, l_error++);
                }
            }
            cerr << "...done" << endl;
        }
        catch (const exception& e) {
            cerr << e.what() << endl;
        }
    }
}

void ResultsCalculator::integrate_hard_factor(const Context& ctx, const ThreadLocalContext& tlctx, HardFactorList& hflist, double* l_real, double* l_imag, double* l_error) {
    Integrator* integrator = new Integrator(&ctx, &tlctx, hflist);
    if (trace) {
        integrator->set_callback(write_data_point);
    }
    else if (minmax) {
        integrator->set_callback(store_minmax);
    }
    integrator->set_cubature_callback(cubature_eprint_callback);
    integrator->set_miser_callback(miser_eprint_callback);
    integrator->set_vegas_callback(vegas_eprint_callback);
    integrator->set_quasi_callback(quasi_eprint_callback);
    integrator->integrate(l_real, l_imag, l_error);
    _valid++;
    delete integrator;
}

/**
 * Write the list of results in a ResultsCalculator to the given output stream.
 */
ostream& operator<<(ostream& out, ResultsCalculator& rc) {
    bool all_valid = true;
    _Setw lw = setw(6), rw = setw(26); // "label width" and "result width"
    // write headers
    if (rc.separate) {
        out << lw << left << "pT" << lw << "Y";
        for (size_t hfgindex = 0; hfgindex < rc._hfglen; hfgindex++) {
            out << rw << rc.hfgnames[hfgindex];
            size_t hflen = rc.hfgroups[hfgindex]->size();
            for (size_t hfindex = 1; hfindex < hflen; hfindex++) {
                out << rw << " ";
            }
        }
        out << rw << "total" << endl;
        out << lw << " " << lw << " ";
        for (vector<string>::iterator termname_iterator = rc.hfnames.begin(); termname_iterator != rc.hfnames.end(); termname_iterator++) {
            out << rw << *termname_iterator;
        }
        out << endl;
    }
    else {
        out << lw << left << "pT" << lw << "Y";
        for (vector<string>::iterator it = rc.hfgnames.begin(); it != rc.hfgnames.end(); it++) {
            out << rw << *it;
        }
        out << rw << "total" << endl;
    }
    
    // write data
    double l_real, l_imag, l_error;
    for (size_t ccindex = 0; ccindex < rc.cc.size(); ccindex++) {
        out << lw << sqrt(rc.cc[ccindex].pT2);
        out << lw << rc.cc[ccindex].Y;

        double total = 0;
        size_t hfglen = rc.separate ? rc._hflen : rc._hfglen;
        for (size_t hfgindex = 0; hfgindex < hfglen; hfgindex++) {
            if (rc.valid(ccindex, hfgindex)) {
                rc.result(ccindex, hfgindex, &l_real, &l_imag, &l_error);
                ostringstream s;
                s << l_real << "±" << l_error;
                /* needs an extra space because the "±" counts as two chars for computing the
                 * field width, but only displays as one character */
                out << rw << s.str() << " ";
                total += l_real;
            }
            else {
                out << rw << "---";
                all_valid = false;
            }
        }
        if (all_valid) {
            out << rw << total << endl;
        }
        else {
            out << rw << "---" << endl;
        }
    }
    if (!all_valid) {
        out << "WARNING: some results were not computed" << endl;
    }

    if (rc.minmax) {
        out << "xx\t" << min_ictx.xx << "\t" << max_ictx.xx << "\t" << endl;
        out << "xy\t" << min_ictx.xy << "\t" << max_ictx.xy << "\t" << endl;
        out << "yx\t" << min_ictx.yx << "\t" << max_ictx.yx << "\t" << endl;
        out << "yy\t" << min_ictx.yy << "\t" << max_ictx.yy << "\t" << endl;
        out << "bx\t" << min_ictx.bx << "\t" << max_ictx.bx << "\t" << endl;
        out << "by\t" << min_ictx.by << "\t" << max_ictx.by << "\t" << endl;
        out << "q1x\t" << min_ictx.q1x << "\t" << max_ictx.q1x << "\t" << endl;
        out << "q1y\t" << min_ictx.q1y << "\t" << max_ictx.q1y << "\t" << endl;
        out << "q2x\t" << min_ictx.q2x << "\t" << max_ictx.q2x << "\t" << endl;
        out << "q2y\t" << min_ictx.q2y << "\t" << max_ictx.q2y << "\t" << endl;
        out << "q3x\t" << min_ictx.q3x << "\t" << max_ictx.q3x << "\t" << endl;
        out << "q3y\t" << min_ictx.q3y << "\t" << max_ictx.q3y << "\t" << endl;
        out << "z\t" << min_ictx.z << "\t" << max_ictx.z << "\t" << endl;
        out << "xi\t" << min_ictx.xi << "\t" << max_ictx.xi << "\t" << endl;
        out << "xip\t" << min_ictx.xiprime << "\t" << max_ictx.xiprime << "\t" << endl;
        out << "kT\t" << min_ictx.kT << "\t" << max_ictx.kT << "\t" << endl;
        out << "kT2\t" << min_ictx.kT2 << "\t" << max_ictx.kT2 << "\t" << endl;
        out << "xp\t" << min_ictx.xp << "\t" << max_ictx.xp << "\t" << endl;
        out << "xg\t" << min_ictx.xg << "\t" << max_ictx.xg << "\t" << endl;
        out << "Qs2\t" << min_ictx.Qs2 << "\t" << max_ictx.Qs2 << "\t" << endl;
        out << "alphas\t" << min_ictx.alphas << "\t" << max_ictx.alphas << "\t" << endl;
    }
}

/**
 * An exception to be thrown when there is an error in the GSL code.
 */
class GSLException : public exception {
private:
    string _reason;
    string _file;
    int _line;
    int _gsl_errno;
    string _message;
public:
    GSLException(const char* reason, const char* file, int line, int gsl_errno) throw() :
        _reason(reason), _file(file), _line(line), _gsl_errno(gsl_errno) {
        ostringstream s;
        s << "GSL error " << gsl_errno << "(" << gsl_strerror(gsl_errno) << "): " << reason << " at " << file << ":" << line;
        _message = s.str();
    }
    GSLException(const GSLException& e) throw() :
        _reason(e._reason), _file(e._file), _line(e._line), _gsl_errno(e._gsl_errno), _message(e._message) {
    }
    GSLException& operator=(const GSLException& e) throw() {
        _reason = e._reason;
        _file = e._file;
        _line = e._line;
        _gsl_errno = e._gsl_errno;
        _message = e._message;
    }
    ~GSLException() throw() {
    }
    const string& reason() const {
        return _reason;
    }
    const string& file() const {
        return _file;
    }
    const int line() const {
        return _line;
    }
    const int gsl_errno() const {
        return _gsl_errno;
    }
    const char* what() const throw() {
        return _message.c_str();
    }
};

/**
 * GSL error handler function that throws a GSLException.
 */
void gsl_error_throw(const char* reason, const char* file, int line, int gsl_errno) {
    throw GSLException(reason, file, line, gsl_errno);
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
    _exit(2);
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
    cc.create_contexts();
    
    ThreadLocalContext tlctx(cc);

    /* First write out all the configuration variables. Having the configuration written
     * out as part of the output file makes it easy to tell what parameters were used in
     * and given run, and is also useful in case we want to reproduce a run.
     */
#ifdef GIT_REVISION
    cout << "# git revision " << GIT_REVISION << endl;
#endif
    cout << cc << "------------" << endl;

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
    catch (const exception& e) {
        cerr << e.what() << endl;
        return 1;
    }
    catch (const char* c) {
        cerr << c << endl;
        return 1;
    }
    catch (...) {
        cerr << "unknown error" << endl;
        return 1;
    }
}
