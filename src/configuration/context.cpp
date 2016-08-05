/*
 * Part of SOLO
 *
 * Copyright 2012-2015 David Zaslavsky
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
#include <cstdlib>
#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <gsl/gsl_math.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>
#include "context.h"
#include "../gluondist/gluondist.h"
#include "../integration/integrator.h"
#include "../log.h"
#include "../utils/utils.h"

using namespace std;

#define check_property_default(p, typename, parse, default) \
    itit = m_config.equal_range(canonicalize(#p));\
    typename p;\
    if (itit.first == itit.second) {\
        p = default;\
        m_config.set(#p, format(p));\
    }\
    else {\
        p = parse(itit);\
    }
#define check_property(p, typename, parse) \
    itit = m_config.equal_range(canonicalize(#p));\
    typename p;\
    if (itit.first == itit.second) {\
        throw MissingPropertyException(#p);\
    }\
    else {\
        p = parse(itit);\
    }

const string trim_lower(const string& i_str) {
    string str = trim(i_str, " \t\n");
    transform(str.begin(), str.end(), str.begin(), ::tolower);
    return str;
}

string canonicalize(const string& i_key) {
    const string key = trim_lower(i_key);
    if (key == "lambda_qcd") {
        return "lambdaqcd";
    }
    else if (key == "lambda_mv") {
        return "lambdamv";
    }
    else if (key == "mu^2") {
        return "mu2";
    }
    else if (key == "ncolors") {
        return "nc";
    }
    else if (key == "nflavors") {
        return "nf";
    }
    else if (key == "colorfactor") {
        return "cf";
    }
    else if (key == "hf_definitions") {
        return "hardfactor_definitions";
    }
    else if (key == "alpha_s") {
        return "alphas";
    }
    else if (key == "coupling" || key == "cpl") {
        return "coupling_type";
    }
    else if (key == "gdist" || key == "gluon_distribution" || key == "gluon distribution" || key == "gluon dist") {
        return "gdist_type";
    }
    else if (key == "quasirandom generator type" || key == "qrng type") {
        return "quasirandom_generator_type";
    }
    else if (key == "pseudorandom generator type" || key == "rng type") {
        return "pseudorandom_generator_type";
    }
    else if (key == "pseudorandom generator seed" || key == "rng seed" || key == "seed") {
        return "pseudorandom_generator_seed";
    }
    else if (key == "qmc iterations" || key == "quasi iterations") {
        return "quasi_iterations";
    }
    else if (key == "qmc absolute error" || key == "quasi absolute error" || key == "absolute error") {
        return "abserr";
    }
    else if (key == "qmc relative error" || key == "quasi relative error" || key == "relative error") {
        return "relerr";
    }
    else if (key == "integration strategy" || key == "strategy") {
        return "integration_strategy";
    }
    else if (key == "c" || key == "centrality parameter") {
        return "centrality";
    }
    else if (key == "resum_c" || key == "resummation constant") {
        return "resummation_constant";
    }
    else if (key == "a") {
        return "mass_number";
    }
    else if (key == "exact kinematics" || key == "exact") {
        return "exact_kinematics";
    }
    else {
        return key;
    }
}

template<typename T>
const string format(const T& v) {
    ostringstream s;
    s << v;
    return s.str();
}

static inline void assert_range_size_one(const pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    multimap<string, string>::iterator el1 = range.first;
    multimap<string, string>::iterator el2 = range.second;
    assert(++el1 == el2);
}

size_t parse_size(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    assert_range_size_one(range);
    return static_cast<size_t>(strtoul(range.first->second.c_str(), NULL, 0));
}

unsigned long int parse_ulong(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    assert_range_size_one(range);
    return strtoul(range.first->second.c_str(), NULL, 0);
}

double parse_double(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    assert_range_size_one(range);
    return strtod(range.first->second.c_str(), NULL);
}

bool parse_boolean(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    assert_range_size_one(range);
    return range.first->second == "true"
        || range.first->second == "yes"
        || range.first->second == "1";
}

string& parse_string(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    assert_range_size_one(range);
    return range.first->second;
}

projectile_type parse_projectile_type(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    assert_range_size_one(range);
    string val = range.first->second;
    transform(val.begin(), val.end(), val.begin(), ::tolower);
    if (val == "proton") {
        return proton;
    }
    else if (val == "deuteron") {
        return deuteron;
    }
    else {
        GSL_ERROR_VAL("unknown method", GSL_EINVAL, proton);
    }
}

DSSpiNLO::hadron parse_hadron(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    assert_range_size_one(range);
    string val = trim_lower(range.first->second);
    if (val == "pim" || val == "pi-" || val == "pi minus") {
        return DSSpiNLO::pi_minus;
    }
    else if (val == "piz" || val == "pi0" || val == "pi zero") {
        return DSSpiNLO::pi_zero;
    }
    else if (val == "pip" || val == "pi+" || val == "pi plus") {
        return DSSpiNLO::pi_plus;
    }
    else {
        GSL_ERROR_VAL("unknown hadron", GSL_EINVAL, DSSpiNLO::pi_minus);
    }
}

integration_strategy parse_strategy(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    assert_range_size_one(range);
    string val = trim_lower(range.first->second);
    if (val == "plain") {
        return MC_PLAIN;
    }
    else if (val == "miser") {
        return MC_MISER;
    }
    else if (val == "vegas") {
        return MC_VEGAS;
    }
    else if (val == "quasi") {
        return MC_QUASI;
    }
    else {
        GSL_ERROR_VAL("unknown method", GSL_EINVAL, MC_VEGAS);
    }
}

const gsl_rng_type* parse_rng_type(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    assert_range_size_one(range);
    static const gsl_rng_type** types = gsl_rng_types_setup();
    for (const gsl_rng_type** t = types; *t != NULL; t++) {
        string name((*t)->name);
        if (name == range.first->second) {
            return *t;
        }
    }
    GSL_ERROR_VAL("unknown generator", GSL_EINVAL, NULL);
}

const gsl_qrng_type* parse_qrng_type(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    assert_range_size_one(range);
    string val = trim_lower(range.first->second);
    if (val == "niederreiter_2") {
        return gsl_qrng_niederreiter_2;
    }
    else if (val == "sobol") {
        return gsl_qrng_sobol;
    }
    else if (val == "halton") {
        return gsl_qrng_halton;
    }
    else if (val == "reversehalton") {
        return gsl_qrng_reversehalton;
    }
    else {
        GSL_ERROR_VAL("unknown generator", GSL_EINVAL, NULL);
    }
}

vector<double> parse_double_vector(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    vector<double> v;
    for (multimap<string, string>::iterator it = range.first; it != range.second; it++) {
        v.push_back(strtod(it->second.c_str(), NULL));
    }
    return v;
}

vector<unsigned long int> parse_ulong_vector(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    vector<unsigned long int> v;
    for (multimap<string, string>::iterator it = range.first; it != range.second; it++) {
        v.push_back(strtoul(it->second.c_str(), NULL, 0));
    }
    return v;
}

vector<string> parse_string_vector(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    vector<string> v;
    for (multimap<string, string>::iterator it = range.first; it != range.second; it++) {
        v.push_back(it->second);
    }
    return v;
}

struct Mu2Strategy {
    double _mu2;
    Mu2Strategy(double mu2) : _mu2(mu2) {}
    virtual double mu2(double pT2) {
        return _mu2;
    }
};
struct PT2x4Mu2Strategy {
    virtual double mu2(double pT2) {
        return 4*pT2;
    }
};


ContextCollection::ContextCollection(const Configuration& conf) :
  m_config(conf),
  m_gdist(NULL),
  m_cpl(NULL),
  m_fs(NULL),
  trace_gdist(false) {
    create_contexts();
}

ContextCollection::~ContextCollection() {
    delete m_gdist;
    m_gdist = NULL;
    delete m_cpl;
    m_cpl = NULL;
}

GBWGluonDistribution* ContextCollection::create_gbw_gluon_distribution() {
    return new GBWGluonDistribution(Q02, x0, lambda);
}

MVGluonDistribution* ContextCollection::create_mv_gluon_distribution() {
    pair<multimap<string, string>::iterator, multimap<string, string>::iterator> itit;
    double Ymin = min(Y);
    double Ymax = max(Y);
    double pTmin = min(pT);
    check_property_default(lambdaMV, double, parse_double, 0.241)
    check_property_default(gammaMV,  double, parse_double, 1)
    check_property_default(q2minMV,  double, parse_double, 1e-6)
    /* q2maxMV by default is the lesser of
     *  (2 qxmax + sqrt(smax) / exp(Ymin))^2 + (2 qymax)^2
     * or 1e4. The value is empirically chosen to be small enough not to
     * produce numerical errors in the MV distribution's Fourier transform,
     * but large enough that the distribution above q2maxMV follows
     * 1/q^4 scaling. This really should be tuned in the configuration file,
     * as the ideal value may differ from one calculation to the next.
     */
    check_property_default(q2maxMV,  double, parse_double, min(gsl_pow_2(2 * inf + sqs / exp(Ymin)) + gsl_pow_2(2 * inf), 1.e4))
    check_property_default(YminMV, double, parse_double, 2 * Ymin)
    check_property_default(YmaxMV, double, parse_double, Ymax - log(pTmin) + log(sqs))
    check_property_default(gdist_subinterval_limit, size_t, parse_size, 10000)
    logger << "Creating MV gluon distribution with " << q2minMV << " < k2 < " << q2maxMV << ", " << YminMV << " < Y < " << YmaxMV << endl;
    return new MVGluonDistribution(lambdaMV, gammaMV, q2minMV, q2maxMV, YminMV, YmaxMV, Q02, x0, lambda, gdist_subinterval_limit);
}

FixedSaturationMVGluonDistribution* ContextCollection::create_fmv_gluon_distribution() {
    pair<multimap<string, string>::iterator, multimap<string, string>::iterator> itit;
    double Ymin = min(Y);
    check_property_default(lambdaMV, double, parse_double, 0.241)
    check_property_default(gammaMV,  double, parse_double, 1)
    check_property_default(q2minMV,  double, parse_double, 1e-6)
    // q2max = (2 qxmax + sqrt(smax) / exp(Ymin))^2 + (2 qymax)^2
    check_property_default(q2maxMV,  double, parse_double, gsl_pow_2(2 * inf + sqs / exp(Ymin)) + gsl_pow_2(2 * inf))
    check_property(YMV, double, parse_double)
    logger << "Creating fMV gluon distribution with " << q2minMV << " < k2 < " << q2maxMV << ", Y = " << YMV << endl;
    check_property_default(gdist_subinterval_limit, size_t, parse_size, 10000)
    return new FixedSaturationMVGluonDistribution(lambdaMV, gammaMV, q2minMV, q2maxMV, YMV, Q02, x0, lambda, gdist_subinterval_limit);
}

PlateauPowerGluonDistribution* ContextCollection::create_pp_gluon_distribution() {
    pair<multimap<string, string>::iterator, multimap<string, string>::iterator> itit;
    double Ymin = min(Y);
    double Ymax = max(Y);
    double pTmin = min(pT);
    check_property_default(gammaPP, double, parse_double, 4)
    check_property_default(r2minPP, double, parse_double, 1e-6)
    check_property_default(r2maxPP, double, parse_double, 2 * gsl_pow_2(2 * inf))
    check_property_default(YminPP,  double, parse_double, 2 * Ymin)
    check_property_default(YmaxPP,  double, parse_double, Ymax - log(pTmin) + log(sqs))
    check_property_default(gdist_subinterval_limit, size_t, parse_size, 10000)
    logger << "Creating plateau-power gluon distribution with " << r2minPP << " < r2 < " << r2maxPP << ", " << YminPP << " < Y < " << YmaxPP << endl;
    return new PlateauPowerGluonDistribution(gammaPP, r2minPP, r2maxPP, YminPP, YmaxPP, Q02, x0, lambda, gdist_subinterval_limit);
}

FileDataGluonDistribution* ContextCollection::create_file_gluon_distribution(GluonDistribution* lower_dist = NULL, GluonDistribution* upper_dist = NULL, const bool extended = false) {
    pair<multimap<string, string>::iterator, multimap<string, string>::iterator> itit;
    // check that if creating an extended distribution, we're not passing
    // anything for lower_dist or upper_dist, just to make sure there are
    // no programmer misconceptions about what the method is doing
    assert(extended || (lower_dist == NULL && upper_dist == NULL));

    check_property(gdist_position_filename, string, parse_string)
    check_property(gdist_momentum_filename, string, parse_string)
    check_property_default(satscale_source, string, parse_string, "analytic")
    check_property_default(xinit, double, parse_double, 0.01)
    logger << "Reading gluon distribution from " << gdist_position_filename << " (pos) and " << gdist_momentum_filename << " (mom)" << endl;

    if (satscale_source == "analytic") {
        if (extended) {
            return new ExtendedFileDataGluonDistribution(
             gdist_position_filename, gdist_momentum_filename,
             Q02, x0, lambda, xinit,
             lower_dist, upper_dist);
        }
        else {
            return new FileDataGluonDistribution(
             gdist_position_filename, gdist_momentum_filename,
             Q02, x0, lambda, xinit);
        }
    }
    else if (satscale_source == "extract from momentum" || satscale_source == "extract from position") {
        check_property(satscale_threshold, double, parse_double)
        FileDataGluonDistribution::satscale_source_type satscale_source_constant;
        if (satscale_source == "extract from momentum") {
            satscale_source_constant = FileDataGluonDistribution::MOMENTUM_THRESHOLD;
        }
        else if (satscale_source == "extract from position") {
            satscale_source_constant = FileDataGluonDistribution::POSITION_THRESHOLD;
        }
        else {
            assert(false);
        }
        if (extended) {
            return new ExtendedFileDataGluonDistribution(
             gdist_position_filename, gdist_momentum_filename,
             xinit, satscale_source_constant, satscale_threshold,
             lower_dist, upper_dist);
        }
        else {
            return new FileDataGluonDistribution(
             gdist_position_filename, gdist_momentum_filename,
             xinit, satscale_source_constant, satscale_threshold);
        }
    }
    else {
        throw InvalidPropertyValueException<string>("satscale_source", satscale_source);
    }
}

GluonDistribution* ContextCollection::create_gluon_distribution(const string& gdist_type) {
    if (gdist_type == "gbw") {
        return create_gbw_gluon_distribution();
    }
    else if (gdist_type == "mv") {
        return create_mv_gluon_distribution();
    }
    else if (gdist_type == "fmv") {
        return create_fmv_gluon_distribution();
    }
    else if (gdist_type == "plateau-power" || gdist_type == "pp") {
        return create_pp_gluon_distribution();
    }
    else if (gdist_type == "file") {
        return create_file_gluon_distribution();
    }
    else {
        vector<string> pieces = split(gdist_type, "+");
        GluonDistribution* lower_dist = NULL;
        GluonDistribution* upper_dist = NULL;
        if (pieces.size() == 1) {
            // gdist_type = "efile"
            if (pieces[0] != "efile") {
                throw InvalidPropertyValueException<string>("gdist_type", gdist_type);
            }
        }
        if (pieces.size() == 2) {
            // gdist_type = "<type>+file" or "file+<type>"
            if (pieces[0] == "file" || pieces[0] == "efile") {
                if (pieces[1] == "file" || pieces[1] == "efile") {
                    throw InvalidPropertyValueException<string>("gdist_type", gdist_type);
                }
                upper_dist = create_gluon_distribution(pieces[1]);
            }
            else if (pieces[1] == "file" || pieces[1] == "efile") {
                lower_dist = create_gluon_distribution(pieces[0]);
            }
            else {
                throw InvalidPropertyValueException<string>("gdist_type", gdist_type);
            }
        }
        else if (pieces.size() == 3) {
            // gdist_type = "<type>+file+<type>"
            if (pieces[1] == "file" || pieces[1] == "efile") {
                if (pieces[0] == "file" || pieces[0] == "efile" || pieces[2] == "file" || pieces[2] == "efile") {
                    throw InvalidPropertyValueException<string>("gdist_type", gdist_type);
                }
                lower_dist = create_gluon_distribution(pieces[0]);
                upper_dist = create_gluon_distribution(pieces[2]);
            }
            else {
                throw InvalidPropertyValueException<string>("gdist_type", gdist_type);
            }
        }
        else {
            throw InvalidPropertyValueException<string>("gdist_type", gdist_type);
        }
        return create_file_gluon_distribution(lower_dist, upper_dist, true);
    }
}

void ContextCollection::create_contexts() {
    bool _c0r_optimization = false;
    pair<multimap<string, string>::iterator, multimap<string, string>::iterator> itit;

    check_property_default( x0,           double, parse_double, 0.000304)
    check_property(         mass_number,  double, parse_double)
    check_property(         centrality,   double, parse_double)
    check_property_default( lambda,       double, parse_double, 0.288)
    check_property_default( Nc,           double, parse_double, 3)
    check_property_default( Nf,           double, parse_double, 3)
    check_property_default( CF,           double, parse_double, 1.5)
    check_property_default( TR,           double, parse_double, 0.5)
    check_property(         Sperp,        double, parse_double)
    check_property(         sqs,          double, parse_double)
    check_property(         pT,           vector<double>, parse_double_vector)
    check_property(         Y,            vector<double>, parse_double_vector)
    check_property(         projectile,   projectile_type, parse_projectile_type)
    check_property(         hadron,       DSSpiNLO::hadron, parse_hadron)
    check_property(         hardfactor_definitions, vector<string>, parse_string_vector)
    check_property_default( pdf_filename, string, parse_string, "mstw2008nlo.00.dat")
    check_property_default( ff_filename,  string, parse_string, "PINLO.DAT")
    check_property_default( integration_strategy, integration_strategy, parse_strategy, MC_VEGAS)
    check_property_default( cubature_iterations, size_t, parse_size, 1000000)
    check_property_default( miser_iterations, size_t, parse_size, 1000000)
    check_property_default( vegas_initial_iterations, size_t, parse_size, 100000)
    check_property_default( vegas_incremental_iterations, size_t, parse_size, 100000)
    check_property_default( quasi_iterations, size_t, parse_size, 1000000)
    check_property_default( abserr, double, parse_double, 1e-20)
    check_property_default( relerr, double, parse_double, 0)
    check_property_default( css_r_regularization, bool, parse_boolean, false)
    check_property_default( css_r2_max, double, parse_double, 0) // 0 is a dummy value for when the regularization is not being used
    check_property_default( resummation_constant, double, parse_double, 1)
    check_property_default( xif, double, parse_double, 0)
    check_property_default( exact_kinematics, bool, parse_boolean, false)

    // Adapted from the GSL source code - basically this reimplements gsl_rng_env_setup
    if (!m_config.contains("quasirandom_generator_type")) {
        const char* qtype = getenv("GSL_QRNG_TYPE");
        m_config.set("quasirandom_generator_type", qtype == NULL ? "halton" : qtype);
    }
    if (!m_config.contains("pseudorandom_generator_type")) {
        const char* type = getenv("GSL_RNG_TYPE");
        m_config.set("pseudorandom_generator_type", type == NULL ? "mt19937" : type);
    }
    if (!m_config.contains("pseudorandom_generator_seed")) {
        const char* seed = getenv("GSL_RNG_SEED");
        m_config.set("pseudorandom_generator_seed", seed == NULL ? "0" : seed);
    }

    check_property(         quasirandom_generator_type, const gsl_qrng_type*, parse_qrng_type)
    check_property(         pseudorandom_generator_type, const gsl_rng_type*, parse_rng_type)
    check_property(         pseudorandom_generator_seed, vector<unsigned long int>, parse_ulong_vector)

    /*
     * For dimensions in which the lower and upper bounds are -infinity and +infinity,
     * we have to pick a finite value to cut off the integral. It shouldn't be made
     * too large, though, because if it is, the Monte Carlo sampling will miss the
     * peak in the integrand entirely and just output zero all the time.
     */
    check_property_default( inf, double, parse_double, 40)

    /* The default cutoff is zero, since most integrands are well-behaved at 0. */
    check_property_default( cutoff, double, parse_double, 0)

    this->Q02 = centrality * pow(mass_number, 1./3.);
    this->x0 = x0;
    this->lambda = lambda;
    this->pT = pT;
    this->Y = Y;
    this->sqs = sqs;
    this->inf = inf;

    if (css_r_regularization && css_r2_max <= 0) {
        GSL_ERROR_VOID("invalid r_max", GSL_EINVAL);
    }

    // create gluon distribution
    assert (m_gdist == NULL);
    check_property(gdist_type, string, parse_string)
    gdist_type = trim_lower(gdist_type);
    m_gdist = create_gluon_distribution(gdist_type);
    assert(m_gdist != NULL);
    if (trace_gdist) {
        logger << "Activating gluon distribution trace wrapper" << endl;
        m_gdist = new GluonDistributionTraceWrapper(m_gdist);
    }
    assert(m_gdist != NULL);

    // create coupling
    assert(m_cpl == NULL);
    check_property(coupling_type, string, parse_string)
    if (coupling_type == "fixed") {
        check_property_default(alphas, double, parse_double, 0.2)
        m_cpl = new FixedCoupling(alphas);
    }
    else if (coupling_type == "running") {
        check_property_default(lambdaQCD, double, parse_double, sqrt(0.0588))
        check_property_default(regulator, double, parse_double, 1.0)
        check_property_default(Ncbeta, double, parse_double, (11.0 * Nc - 2.0 * Nf) / 12.0)
        m_cpl = new LORunningCoupling(lambdaQCD, Ncbeta, regulator);
    }
    else {
        throw InvalidPropertyValueException<string>("coupling_type", coupling_type);
    }
    assert(m_cpl != NULL);

    // create factorization scale strategy
    assert(m_fs == NULL);
    check_property_default(factorization_scale, string, parse_string, "fixed")
    factorization_scale = trim_lower(factorization_scale);
    if (factorization_scale == "fixed") {
        /* special case: check for old config file format with
         *  mu2 = 4pT2
         * and convert it to new format
         *  factorization_scale = 4pT2
         */
        itit = m_config.equal_range(canonicalize("mu2"));
        if (trim_lower(itit.first->second) == "4pt2") {
            m_fs = new PTProportionalFactorizationScale(4);
        }
        else {
            // this is the normal case
            check_property_default(mu2, double, parse_double, 10)
            m_fs = new FixedFactorizationScale(mu2);
        }
    }
    else if (factorization_scale == "4pt2") {
        m_fs = new PTProportionalFactorizationScale(4);
    }
    else if (factorization_scale == "cpt2") {
        check_property(factorization_scale_coefficient, double, parse_double)
        m_fs = new PTProportionalFactorizationScale(factorization_scale_coefficient);
    }
    else if (factorization_scale == "c0r") {
        check_property_default(c0r_optimization, bool, parse_boolean, true)
        _c0r_optimization = c0r_optimization;
        m_fs = new RPerpFactorizationScale(4 * exp(-2*M_EULER)); // this value is c_0^2, with c_0 defined in the long paper
    }
    assert(m_fs != NULL);

    // create contexts
    for (vector<double>::iterator pTit = pT.begin(); pTit != pT.end(); pTit++) {
        for (vector<double>::iterator Yit = Y.begin(); Yit != Y.end(); Yit++) {
            for (vector<unsigned long int>::iterator seedit = pseudorandom_generator_seed.begin(); seedit != pseudorandom_generator_seed.end(); seedit++) {
                double pT = *pTit;
                double Y = *Yit;
                unsigned long int seed = *seedit;
                try {
                    const Context c = {
                      x0, mass_number, centrality, lambda, Nc, Nf, CF, TR, Sperp,
                      gsl_pow_2(pT), sqs, Y,
                      hardfactor_definitions, pdf_filename, ff_filename,
                      quasirandom_generator_type, pseudorandom_generator_type, seed,
                      m_gdist, m_cpl, m_fs, _c0r_optimization, css_r_regularization, css_r2_max,
                      resummation_constant,
                      xif,
                      exact_kinematics,
                      projectile, hadron,
                      integration_strategy,
                      abserr, relerr,
                      cubature_iterations, miser_iterations,
                      vegas_initial_iterations, vegas_incremental_iterations, quasi_iterations,
                      inf, cutoff,
                      Context::compute_Q02x0lambda(centrality, mass_number, x0, lambda),
                      Context::compute_tau(pT, sqs, Y),
                      Context::compute_tauhat(pT, sqs, Y)
                    };
                    c.check_kinematics();
                    push_back(c);
                }
                catch (const InvalidKinematicsException& e) {
                    logger << "Failed to create context at pT = " << pT << ", Y = " << Y << ": " << e.what() << endl;
                }
            }
        }
    }
}

const Configuration& ContextCollection::config() const {
    return m_config;
}

std::ostream& operator<<(std::ostream& out, const integration_strategy& strategy) {
    switch (strategy) {
        case MC_PLAIN:
            out << "plain";
            break;
        case MC_MISER:
            out << "miser";
            break;
        case MC_VEGAS:
            out << "vegas";
            break;
        case MC_QUASI:
            out << "quasi";
            break;
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, const projectile_type& proj) {
    switch (proj) {
        case proton:
            out << "proton";
            break;
        case deuteron:
            out << "deuteron";
            break;
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, const DSSpiNLO::hadron& hadron) {
    switch (hadron) {
        case DSSpiNLO::pi_minus:
            out << "pi-";
            break;
        case DSSpiNLO::pi_zero:
            out << "pi0";
            break;
        case DSSpiNLO::pi_plus:
            out << "pi+";
            break;
    }
    return out;
}

template<typename T>
ostream& operator<<(ostream& out, const vector<T>& vec) {
    for (typename vector<T>::const_iterator it = vec.begin(); it != vec.end(); it++) {
        if (it != vec.begin()) {
            out << ", ";
        }
        out << *it;
    }
    return out;
}

std::ostream& operator<<(std::ostream& out, const Context& ctx) {
    out << "x0\t= "         << ctx.x0           << endl;
    out << "mass number\t= "<< ctx.mass_number  << endl;
    out << "centrality\t= " << ctx.centrality   << endl;
    out << "lambda\t= "     << ctx.lambda       << endl;
    out << "Nc\t= "         << ctx.Nc           << endl;
    out << "Nf\t= "         << ctx.Nf           << endl;
    out << "CF\t= "         << ctx.CF           << endl;
    out << "TR\t= "         << ctx.TR           << endl;
    out << "Sperp\t= "      << ctx.Sperp        << endl;
    out << "pT2\t= "        << ctx.pT2          << endl;
    out << "sqs\t= "        << ctx.sqs          << endl;
    out << "Y\t= "          << ctx.Y            << endl;
    out << "hf_definitions\t= " << ctx.hardfactor_definitions << endl;
    out << "pdf_filename\t= " << ctx.pdf_filename << endl;
    out << "ff_filename\t= " << ctx.ff_filename  << endl;
    out << "projectile\t= " << ctx.projectile << endl;
    out << "hadron\t= " << ctx.hadron << endl;
    out << "integration_strategy\t= " << ctx.strategy << endl;
    out << "cubature_iterations\t= " << ctx.cubature_iterations << endl;
    out << "miser_iterations\t= " << ctx.miser_iterations << endl;
    out << "vegas_initial_iterations\t= " << ctx.vegas_initial_iterations << endl;
    out << "vegas_incremental_iterations\t= " << ctx.vegas_incremental_iterations << endl;
    out << "quasi_iterations\t= " << ctx.quasi_iterations << endl;
    out << "abserr\t= " << ctx.abserr << endl;
    out << "relerr\t= " << ctx.relerr << endl;
    out << "inf\t= " << ctx.inf << endl;
    out << "cutoff\t= " << ctx.cutoff << endl;
    out << "gluon distribution\t = " << ctx.gdist << endl;
    out << "coupling\t = " << ctx.cpl << endl;
    out << "factorization scale\t = " << ctx.fs << endl;
    out << "c0r optimization\t = " << ctx.c0r_optimization << endl;
    out << "CSS r regularization\t = " << ctx.css_r_regularization << endl;
    out << "CSS r_max\t = " << ctx.css_r2_max << endl;
    out << "resummation constant\t = " << ctx.resummation_constant << endl;
    out << "rapidity factorization scale\t = " << ctx.xif << endl;
    out << "exact kinematics\t = " << ctx.exact_kinematics << endl;
    out << "quasirandom generator type: " <<  ctx.quasirandom_generator_type->name << endl;
    out << "pseudorandom generator type: " <<  ctx.pseudorandom_generator_type->name << endl;
    out << "pseudorandom generator seed: " << ctx.pseudorandom_generator_seed << endl;
    return out;
}

void Context::check_kinematics() const {
    if (compute_Q02x0lambda(centrality, mass_number, x0, lambda) != Q02x0lambda) {
        throw InvalidPropertyValueException<double>("Q02x0lambda", Q02x0lambda, "value provided in initializer does not match the one calculated from other context parameters");
    }
    if (compute_tau(sqrt(pT2), sqs, Y) != tau) {
        throw InvalidPropertyValueException<double>("tau", tau, "value provided in initializer does not match the one calculated from other context parameters");
    }
    if (compute_tauhat(sqrt(pT2), sqs, Y) != tauhat) {
        throw InvalidPropertyValueException<double>("tauhat", tauhat, "value provided in initializer does not match the one calculated from other context parameters");
    }
    if (tau > 1) {
        throw InvalidKinematicsException("τ > 1: empty phase space");
    }
    if (tauhat > 1) {
        throw InvalidKinematicsException("\\hat{τ} > 1: empty phase space");
    }
    if (css_r_regularization) {
        assert(css_r2_max > 0);
    }
}

ThreadLocalContext::ThreadLocalContext(const Context& ctx) :
  pdf_object(new c_mstwpdf(ctx.pdf_filename.c_str())),
  ff_object(new DSSpiNLO(ctx.ff_filename.c_str())) {
}

ThreadLocalContext::ThreadLocalContext(const ContextCollection& cc) :
  pdf_object(NULL), ff_object(NULL) {
    multimap<string,string>::const_iterator it = cc.config().find("pdf_filename");
    if (it == cc.m_config.end()) {
        throw "no PDF filename";
    }
    pdf_object = new c_mstwpdf(it->second.c_str());
    it = cc.m_config.find("ff_filename");
    if (it == cc.m_config.end()) {
        throw "no FF filename";
    }
    ff_object = new DSSpiNLO(it->second.c_str());
}

ThreadLocalContext::~ThreadLocalContext() {
    delete pdf_object;
    delete ff_object;
}
