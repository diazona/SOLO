/*
 * A calculation of the NLO cross section of pA->pion collisions
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
#include <iostream>
#include <sstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include "cubature.h"
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

using namespace std;

const int SUCCESS = 0;

static bool minmax = false;

static integration_strategy strategy = MC_VEGAS;

// callbacks
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

static IntegrationContext min_ictx(NULL, NULL);
static IntegrationContext max_ictx(NULL, NULL);

#define store(property) \
  min_ictx.property = min_ictx.property == 0 ? ictx->property : min(min_ictx.property, ictx->property); \
  max_ictx.property = max_ictx.property == 0 ? ictx->property : max(max_ictx.property, ictx->property);

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
    store(alphasbar);
}

void write_nonzero(const IntegrationContext* ictx, const double real, const double imag) {
    if (real != 0 || imag != 0) {
        cerr << real << "\t" << imag << endl;
    }
}

void vegas_eprint_callback(double* p_result, double* p_abserr, gsl_monte_vegas_state* s) {
    cerr << "VEGAS output: " << *p_result << " err: " << *p_abserr << " chisq:" << gsl_monte_vegas_chisq(s) << endl;
}
void miser_eprint_callback(double* p_result, double* p_abserr, gsl_monte_miser_state* s) {
    cerr << "MISER output: " << *p_result << " err: " << *p_abserr << endl;
}

class ResultsCalculator {
private:
    ContextCollection& cc;
    ThreadLocalContext& tlctx;
    vector<HardFactorList*> hfgroups;
    vector<string> hfgnames;
    size_t _valid;
    double* real;
    double* imag;
    double* error;
    friend ostream& operator<<(ostream&, ResultsCalculator&);
public:
    const bool trace;
    const bool minmax;
    ResultsCalculator(ContextCollection& cc, ThreadLocalContext& tlctx, vector<HardFactorList*>& hfgroups, vector<string>& hfgnames, bool trace = false, bool minmax = false) :
      cc(cc), tlctx(tlctx), hfgroups(hfgroups), hfgnames(hfgnames), _valid(0), trace(trace), minmax(minmax) {
        size_t hflen = hfgroups.size();
        size_t cclen = cc.size();
        real = new double[hflen * cclen];
        imag = new double[hflen * cclen];
        error = new double[hflen * cclen];
    };
    ~ResultsCalculator() {
        delete[] real;
        delete[] imag;
        delete[] error;
    };
    size_t index_from(size_t ccindex, size_t hfindex) {
        return ccindex * hfgroups.size() + hfindex;
    };
    bool valid(size_t ccindex, size_t hfindex) {
        return index_from(ccindex, hfindex) < _valid;
    }
    void result(size_t ccindex, size_t hfindex, double* real, double* imag, double* error);
    void calculate();
};

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
            for (vector<HardFactorList*>::iterator hit = hfgroups.begin(); hit != hfgroups.end(); hit++) {
                Integrator* integrator = new Integrator(&ctx, &tlctx, strategy, **hit);
                if (trace) {
                    integrator->set_callback(write_data_point);
                }
                else if (minmax) {
                    integrator->set_callback(store_minmax);
                }
                integrator->set_miser_callback(miser_eprint_callback);
                integrator->set_vegas_callback(vegas_eprint_callback);
                integrator->integrate(l_real++, l_imag++, l_error++);
                _valid++;
                delete integrator;
            }
            cerr << "...done" << endl;
        }
        catch (const exception& e) {
            cerr << e.what() << endl;
        }
    }
}

ostream& operator<<(ostream& out, ResultsCalculator& rc) {
    bool all_valid = true;
    out << "pT\tY\t";
    for (vector<string>::iterator it = rc.hfgnames.begin(); it != rc.hfgnames.end(); it++) {
        out << *it << "\t";
    }
    out << "total" << endl;
    double l_real, l_imag, l_error;
    for (size_t ccindex = 0; ccindex < rc.cc.size(); ccindex++) {
        out << sqrt(rc.cc[ccindex].pT2) << "\t";
        out << rc.cc[ccindex].Y << "\t";

        double total = 0;
        for (size_t hfgindex = 0; hfgindex < rc.hfgroups.size(); hfgindex++) {
            if (rc.valid(ccindex, hfgindex)) {
                rc.result(ccindex, hfgindex, &l_real, &l_imag, &l_error);
                out << l_real << "±" << l_error << "\t";
                total += l_real;
            }
            else {
                out << "---\t";
                all_valid = false;
            }
        }
        if (!all_valid) {
            out << "---" << endl;
        }
        else {
            out << total << endl;
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
        out << "alphasbar\t" << min_ictx.alphasbar << "\t" << max_ictx.alphasbar << "\t" << endl;
    }
}

class ProgramConfiguration {
public:
    ProgramConfiguration(int argc, char** argv);
    ~ProgramConfiguration();
    ContextCollection& context_collection() {
        return cc;
    }
    vector<HardFactorList*>& hard_factor_groups() {
        return hfgroups;
    }
    vector<string>& hard_factor_names() {
        return hfgnames;
    }
    
private:
    integration_strategy strategy;
    string gdist_type;
    bool trace;
    bool minmax;
    bool separate;
    ContextCollection cc;
    vector<string> pT;
    vector<HardFactorList*> hfgroups;
    vector<string> hfgnames;
    
    void parse_hf_spec(const string& spec);
};

ProgramConfiguration::ProgramConfiguration(int argc, char** argv) {
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a == "--miser") {
            strategy = MC_MISER;
        }
        else if (a == "--vegas") {
            strategy = MC_VEGAS;
        }
        else if (a == "--trace") {
            trace = true;
        }
        else if (a == "--minmax") {
            minmax = true;
        }
        else if (a == "--separate") {
            separate = true;
        }
        else if (a == "MV" || a == "GBW") {
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
}

ProgramConfiguration::~ProgramConfiguration() {
    for (vector<HardFactorList*>::iterator it = hfgroups.begin(); it != hfgroups.end(); it++) {
        delete *it;
    }
}

static const char* default_lo_spec = "p:h02qq,p:h02gg";
static const char* default_nlo_spec = "p:h12qq,p:h14qq,p:h12gg,m:h12qqbar,m:h16gg,p:h112gq,p:h122gq,m:h14gq,p:h112qg,p:h122qg,m:h14qg";

void ProgramConfiguration::parse_hf_spec(const string& spec) {
    vector<string> hfnames;
    if (spec == "lo") {
        hfnames = split(default_lo_spec, ", ");
    }
    else if (spec == "nlo") {
        hfnames = split(default_nlo_spec, ", ");
    }
    else {
        hfnames = split(spec, ", ");
    }
    HardFactorList* hfobjs = new HardFactorList();
    assert(hfobjs != NULL);
    for (vector<string>::iterator it = hfnames.begin(); it != hfnames.end(); it++) {
        string s = *it;
        const HardFactorRegistry* registry = position::registry::get_instance(); // default
        if (s[1] == ':') {
            switch (s[0]) {
                case 'm':
                    registry = momentum::registry::get_instance();
                    break;
                case 'p':
                    registry = position::registry::get_instance();
                    break;
            }
            s = s.substr(2);
        }
        const HardFactor* hf = registry->get_hard_factor(s);
        if (hf == NULL) {
            logger << "No such hard factor " << s << endl;
        }
        else {
            hfobjs->push_back(hf);
        }
    }
    assert(hfobjs->size() > 0);
    hfgroups.push_back(hfobjs);
    hfgnames.push_back(spec);
}

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

void gsl_error_throw(const char* reason, const char* file, int line, int gsl_errno) {
    throw GSLException(reason, file, line, gsl_errno);
}

ostream& logger = cerr;

static ResultsCalculator* p_rc = NULL;

void termination_handler(int signal) {
    static bool terminated = false;
    if (!terminated) {
        terminated = true;
        cout << *p_rc;
    }
    _exit(2);
}

int run(int argc, char** argv) {
    gsl_rng_env_setup();
    gsl_set_error_handler(&gsl_error_throw);

    ProgramConfiguration pc(argc, argv);
    ContextCollection cc = pc.context_collection();
    if (cc.empty()) {
        logger << "No momenta or no rapidities specified!" << endl;
        return 1;
    }
    
    ThreadLocalContext tlctx(cc);
    cout << cc << "------------" << endl;

    ResultsCalculator rc(cc, tlctx, pc.hard_factor_groups(), pc.hard_factor_names());
    p_rc = &rc;
    
    struct sigaction siga;
    siga.sa_handler = termination_handler;
    struct sigaction oldsiga;
    sigaction(SIGTERM, &siga, &oldsiga);
    sigaction(SIGINT, &siga, &oldsiga);
    rc.calculate();
    sigaction(SIGTERM, &oldsiga, NULL);
    sigaction(SIGINT, &oldsiga, NULL);
    cout << rc;
    return 0;
}

int main(int argc, char** argv) {
    try {
        return run(argc, argv);
    }
    catch (const exception& e) {
        cerr << e.what() << endl;
        return 1;
    }
    catch (...) {
        cerr << "unknown error" << endl;
        return 1;
    }
}