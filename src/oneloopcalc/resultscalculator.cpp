#include <iomanip>
#include <ostream>
#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>
#include "../hardfactors/hardfactor.h"
#include "../hardfactors/hardfactor_parser.h"
#include "../integration/integrator.h"
#include "../integration/integrationcontext.h"
#include "quasimontecarlo.h"
#include "resultscalculator.h"
#include "trace.h"

using std::bitset;
using std::ostream;
using std::vector;
using trace_variable::trace_vars;

namespace trace_variable {
    // declaration is in trace.h
    bitset<trace_variable::COUNT> trace_vars;
}

/* The following functions are callback functions to be used with
 * Integrator.set_callback(). These would be invoked every time
 * the Monte Carlo routine evaluates the function.
 */

/**
 * A callback function that prints out a bunch of kinematic variables.
 */
void write_data_point(const IntegrationContext* ictx, const double real, const double imag) {
    static ofstream trace_stream("trace.output");
    if (ictx == NULL) {
        trace_stream << endl;
        return;
    }
    #define process(v) if (trace_vars[static_cast<size_t>(trace_variable::v)]) { trace_stream << ictx->v << "\t"; }
    #include "../integration/ictx_var_list.inc"
    #undef process
    trace_stream << real << "\t" << imag << "\t";
    trace_stream << endl;
}

#define process(property) double property;

/** A version of IntegrationContext without the methods */
struct IntegrationContextData {
    #include "../integration/ictx_var_list.inc"
};

static IntegrationContextData min_ictx;
static IntegrationContextData max_ictx;

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

ResultsCalculator::ResultsCalculator(const ProgramConfiguration& pc) :
    cc(pc.config()),
    tlctx(cc),
    result_array_len(cc.size()),
    _hfglen(0),
    _hflen(0),
    trace(pc.trace()),
    minmax(pc.minmax()),
    separate(pc.separate()),
    print_integration_progress(pc.print_integration_progress()),
    xg_min(pc.xg_min()),
    xg_max(pc.xg_max())
{
    assert(hfgroups.empty());
    parse_hf_specs(pc.hfspecs());
    _hfglen = hfgroups.size();
    assert(_hfglen > 0);
    if (separate) {
        for (vector<const HardFactorGroup*>::iterator hfgit = hfgroups.begin(); hfgit != hfgroups.end(); hfgit++) {
            _hflen += (*hfgit)->objects.size();
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
void ResultsCalculator::parse_hf_specs(const vector<string>& hfspecs) {
    // parse the hard factor definition files
    HardFactorParser parser(registry);
    for (vector<string>::const_iterator it = cc[0].hardfactor_definitions.begin(); it != cc[0].hardfactor_definitions.end(); it++) {
        parser.parse_file(*it);
    }
    parser.flush_groups();

    // parse hard factor specifications given on the command line
    for (vector<string>::const_iterator it = hfspecs.begin(); it != hfspecs.end(); it++) {
        const HardFactorGroup* hfg;
        const string& spec = *it;
        if (spec.find(":") != string::npos) {
            // includes a colon, so it is a complete hard factor group specification
            hfg = parser.parse_hard_factor_group(spec);
            if (hfg != NULL) {
                registry.add_hard_factor_group(hfg, true);
            }
        }
        else {
            // no colon, so it references either a group specification defined in the file
            // or earlier on the command line
            hfg = registry.get_hard_factor_group(spec);
            if (hfg == NULL) {
                // or an isolated hard factor
                hfg = parser.parse_hard_factor_group(spec);
                if (hfg == NULL) {
                    throw InvalidHardFactorSpecException(spec, "hard factor group not found");
                }
                else {
                    registry.add_hard_factor_group(hfg, true);
                }
            }
        }
        hfgroups.push_back(hfg);
        hfnames.insert(hfnames.end(), hfg->specifications.begin(), hfg->specifications.end());
    }
    assert(!hfgroups.empty());
    assert(hfnames.size() >= hfgroups.size());
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
    for (ContextCollection::const_iterator it = cc.begin(); it != cc.end(); it++) {
        const Context& ctx = *it;
        cerr << "Beginning calculation at pT = " << sqrt(ctx.pT2) << ", Y = " << ctx.Y << endl;
        try {
            hf_index = 0;
            // recall the definition
            // typedef HardFactorList std::vector<const HardFactor*>
            for (vector<const HardFactorGroup*>::iterator hgit = hfgroups.begin(); hgit != hfgroups.end(); hgit++) {
                if (separate) {
                    // go through the hard factors in each group one at a time
                    HardFactorList one_hf;
                    for (HardFactorList::const_iterator hfit = (*hgit)->objects.begin(); hfit != (*hgit)->objects.end(); hfit++) {
                        one_hf.assign(1, *hfit);
                        integrate_hard_factor(ctx, tlctx, one_hf, index_from(cc_index, hf_index));
                        hf_index++;
                    }
                }
                else {
                    integrate_hard_factor(ctx, tlctx, (*hgit)->objects, index_from(cc_index, hf_index));
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
    Integrator integrator(ctx, tlctx, hflist, xg_min, xg_max);
    if (trace) {
        integrator.set_callback(write_data_point);
    }
    else if (minmax) {
        integrator.set_callback(store_minmax);
    }
    if (print_integration_progress) {
        integrator.set_cubature_callback(cubature_eprint_callback);
        integrator.set_miser_callback(miser_eprint_callback);
        integrator.set_vegas_callback(vegas_eprint_callback);
        integrator.set_quasi_callback(quasi_eprint_callback);
    }
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
    int lw = 6;
    int rw = 14; // "label width" and "result width"

    bool multiseed_mode = false;
    {
        unsigned long int first_seed = rc.cc[0].pseudorandom_generator_seed;
        for (size_t i = 1; i < rc.cc.size(); i++) {
            if (rc.cc[i].pseudorandom_generator_seed != first_seed) {
                multiseed_mode = true;
                break;
            }
        }
    }

    // write headers
    if (rc.separate) {
        out << setw(lw) << left << "pT" << OFS << setw(lw) << "Y" << OFS;
        if (multiseed_mode) {
            out << setw(lw) << "seed" << OFS;
        }
        for (size_t hfgindex = 0; hfgindex < rc._hfglen; hfgindex++) {
            out << setw(rw) << rc.hfgroups[hfgindex]->label << OFS;
            size_t hflen = rc.hfgroups[hfgindex]->objects.size();
            for (size_t hfindex = 1; hfindex < 2 * hflen; hfindex++) {
                out << setw(rw) << BLANK << OFS;
            }
        }
        out << setw(rw) << "total" << endl;

        out << setw(lw) << BLANK << OFS << setw(lw) << BLANK << OFS;
        if (multiseed_mode) {
            out << setw(lw) << "seed" << OFS;
        }
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
        if (multiseed_mode) {
            out << setw(lw) << "seed" << OFS;
        }
        for (vector<const HardFactorGroup*>::iterator it = rc.hfgroups.begin(); it != rc.hfgroups.end(); it++) {
            ostringstream valstream;
            valstream << (*it)->label << "-val";
            out << setw(rw) << valstream.str() << OFS;
            ostringstream errstream;
            errstream << (*it)->label << "-err";
            out << setw(rw) << errstream.str() << OFS;
        }
        out << setw(rw) << "total" << endl;
    }

    // write data
    double l_real, l_imag, l_error;
    size_t hfglen = rc.separate ? rc._hflen : rc._hfglen;
    double* counts = NULL;
    double* means  = NULL;
    double* errors = NULL;
    if (multiseed_mode) {
        counts = new double[hfglen];
        means  = new double[hfglen];
        errors = new double[hfglen];
    }
    bool all_valid = true;
    double last_pt = NAN, last_Y = NAN;
    for (size_t ccindex = 0; ccindex < rc.cc.size(); ccindex++) {
        if (multiseed_mode && (last_pt != rc.cc[ccindex].pT2 || last_Y != rc.cc[ccindex].Y)) {
            if (ccindex > 0) {
                out << setw(lw) << "mean" << OFS << setw(lw) << BLANK << OFS << setw(lw) << BLANK << OFS;
                for (size_t hfgindex = 0; hfgindex < hfglen; hfgindex++) {
                    out << setw(rw) << means[hfgindex] << OFS << setw(rw) << BLANK << OFS;
                }
                out << endl;
                out << setw(lw) << "stddev" << OFS << setw(lw) << BLANK << OFS << setw(lw) << BLANK << OFS;
                for (size_t hfgindex = 0; hfgindex < hfglen; hfgindex++) {
                    out << setw(rw) << sqrt(errors[hfgindex])/counts[hfgindex] << OFS << setw(rw) << BLANK << OFS;
                }
                out << endl;
            }
            memset(counts, 0, hfglen * sizeof(double));
            memset(means,  0, hfglen * sizeof(double));
            memset(errors, 0, hfglen * sizeof(double));

            last_pt = rc.cc[ccindex].pT2;
            last_Y = rc.cc[ccindex].Y;
        }

        out << setw(lw) << sqrt(rc.cc[ccindex].pT2) << OFS;
        out << setw(lw) << rc.cc[ccindex].Y << OFS;
        if (multiseed_mode) {
            out << setw(lw) << rc.cc[ccindex].pseudorandom_generator_seed << OFS;
        }

        double total = 0;
        bool row_valid = true;
        for (size_t hfgindex = 0; hfgindex < hfglen; hfgindex++) {
            if (rc.valid(ccindex, hfgindex)) {
                rc.result(ccindex, hfgindex, &l_real, &l_imag, &l_error);
                out << setw(rw) << l_real << OFS << setw(rw) << l_error << OFS;
                total += l_real;

                if (multiseed_mode) {
                    counts[hfgindex]++;
                    double old_mean = means[hfgindex];
                    means[hfgindex] += (l_real - old_mean) / counts[hfgindex];
                    errors[hfgindex] += (l_real - old_mean) * (l_real - means[hfgindex]);
                }
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
    if (multiseed_mode) {
        out << setw(lw) << "mean" << OFS << setw(lw) << BLANK << OFS << setw(lw) << BLANK << OFS;
        for (size_t hfgindex = 0; hfgindex < hfglen; hfgindex++) {
            out << setw(rw) << means[hfgindex] << OFS << setw(rw) << BLANK << OFS;
        }
        out << endl;
        out << setw(lw) << "stddev" << OFS << setw(lw) << BLANK << OFS << setw(lw) << BLANK << OFS;
        for (size_t hfgindex = 0; hfgindex < hfglen; hfgindex++) {
            out << setw(rw) << sqrt(errors[hfgindex])/counts[hfgindex] << OFS << setw(rw) << BLANK << OFS;
        }
        out << endl;
    }
    if (!all_valid) {
        out << "WARNING: some results were not computed" << endl;
    }

    if (multiseed_mode) {
        delete[] counts;
        delete[] means;
        delete[] errors;
    }

    if (rc.minmax) {
        #define process(v) out << #v << "\t" << min_ictx.v << "\t" << max_ictx.v << "\t" << endl;
        #include "../integration/ictx_var_list.inc"
        #undef process
    }
    return out;
}

