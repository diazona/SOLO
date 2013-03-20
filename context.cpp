/*
 * Part of oneloopcalc
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
#include <cstdlib>
#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_qrng.h>
#include <gsl/gsl_rng.h>
#include "context.h"
#include "gluondist.h"
#include "integrator.h"
#include "log.h"
#include "utils.h"

using namespace std;

#define check_property_default(p, typename, parse, default) \
    itit = options.equal_range(canonicalize(#p));\
    typename p;\
    if (itit.first == itit.second) {\
        p = default;\
        set(#p, format(p));\
    }\
    else {\
        p = parse(itit);\
    }
#define check_property(p, typename, parse) \
    itit = options.equal_range(canonicalize(#p));\
    typename p;\
    if (itit.first == itit.second) {\
        throw MissingPropertyException(#p);\
    }\
    else {\
        p = parse(itit);\
    }

const string canonicalize(const string& i_key) {
    string key = trim(i_key, " \t\n");
    transform(key.begin(), key.end(), key.begin(), ::tolower);
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

size_t parse_size(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    multimap<string, string>::iterator el = range.first;
    assert(++el == range.second);
    return (size_t)strtoul(range.first->second.c_str(), NULL, 0);
}

unsigned long int parse_ulong(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    multimap<string, string>::iterator el = range.first;
    assert(++el == range.second);
    return strtoul(range.first->second.c_str(), NULL, 0);
}

double parse_double(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    multimap<string, string>::iterator el = range.first;
    assert(++el == range.second);
    return strtod(range.first->second.c_str(), NULL);
}

string& parse_string(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    multimap<string, string>::iterator el = range.first;
    assert(++el == range.second);
    return range.first->second;
}

projectile_type parse_projectile_type(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    multimap<string, string>::iterator el = range.first;
    assert(++el == range.second);
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
    multimap<string, string>::iterator el = range.first;
    assert(++el == range.second);
    string val = range.first->second;
    transform(val.begin(), val.end(), val.begin(), ::tolower);
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
    multimap<string, string>::iterator el = range.first;
    assert(++el == range.second);
    string val = range.first->second;
    transform(val.begin(), val.end(), val.begin(), ::tolower);
    if (val == "miser") {
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
    multimap<string, string>::iterator el = range.first;
    assert(++el == range.second);
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
    multimap<string, string>::iterator el = range.first;
    assert(++el == range.second);
    if (range.first->second == "niederreiter_2") {
        return gsl_qrng_niederreiter_2;
    }
    else if (range.first->second == "sobol") {
        return gsl_qrng_sobol;
    }
    else if (range.first->second == "halton") {
        return gsl_qrng_halton;
    }
    else if (range.first->second == "reversehalton") {
        return gsl_qrng_reversehalton;
    }
    else {
        GSL_ERROR_VAL("unknown generator", GSL_EINVAL, NULL);
    }
}

vector<double> parse_vector(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    vector<double> v;
    for (multimap<string, string>::iterator it = range.first; it != range.second; it++) {
        v.push_back(strtod(it->second.c_str(), NULL));
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

void ContextCollection::create_contexts() {
    pair<multimap<string, string>::iterator, multimap<string, string>::iterator> itit;

    check_property(x0,           double, parse_double)
    check_property(A,            double, parse_double)
    check_property(c,            double, parse_double)
    check_property(lambda,       double, parse_double)
    check_property(mu2,          double, parse_double)
    bool mu24pT2 = false;
    // special case: check for mu2=4pT2. This could be done with an object if more complicated schemes are needed.
    itit = options.equal_range(canonicalize("mu2"));
    if (itit.first->second == "4pT2") {
        mu24pT2 = true;
    }
    check_property(Nc,           double, parse_double)
    check_property(Nf,           double, parse_double)
    check_property(CF,           double, parse_double)
    check_property(TR,           double, parse_double)
    check_property(Sperp,        double, parse_double)
    check_property(sqs,          double, parse_double)
    check_property(pT,           vector<double>, parse_vector)
    check_property(Y,            vector<double>, parse_vector)
    check_property(projectile,   projectile_type, parse_projectile_type)
    check_property(hadron,       DSSpiNLO::hadron, parse_hadron)
    check_property(pdf_filename, string, parse_string)
    check_property(ff_filename,  string, parse_string)
    check_property(integration_strategy, integration_strategy, parse_strategy)
    check_property(cubature_iterations, size_t, parse_size)
    check_property(miser_iterations, size_t, parse_size)
    check_property(vegas_initial_iterations, size_t, parse_size)
    check_property(vegas_incremental_iterations, size_t, parse_size)
    check_property(quasi_iterations, size_t, parse_size)
    check_property(abserr, double, parse_double)
    check_property(relerr, double, parse_double)
    check_property(quasirandom_generator_type, const gsl_qrng_type*, parse_qrng_type)
    check_property(pseudorandom_generator_type, const gsl_rng_type*, parse_rng_type)
    check_property(pseudorandom_generator_seed, unsigned long int, parse_ulong)
    
    // TODO: allow different choices of saturation scale
    assert(satscale == NULL);
    satscale = new SaturationScale(c * pow(A, 1./3.), x0, lambda);
    
    // create gluon distribution
    assert (gdist == NULL);
    check_property(gdist_type, string, parse_string)
    if (gdist_type == "GBW") {
        gdist = new GBWGluonDistribution(*satscale);
    }
    else if (gdist_type == "MV") {
        double Ymin = min(Y);
        double Ymax = max(Y);
        double pTmin = min(pT);
        check_property_default(lambdaMV, double, parse_double, 0.241)
        check_property_default(gammaMV,  double, parse_double, 1)
        check_property_default(q2minMV,  double, parse_double, 1e-6)
        // q2max = (2 qxmax + sqrt(smax) / exp(Ymin))^2 + (2 qymax)^2
        check_property_default(q2maxMV,  double, parse_double, gsl_pow_2(2 * inf + sqs / exp(Ymin)) + gsl_pow_2(2 * inf))
        check_property_default(YminMV, double, parse_double, 2 * Ymin)
        check_property_default(YmaxMV, double, parse_double, Ymax - log(pTmin) + log(sqs))
        assert(q2minMV < q2maxMV);
        assert(YminMV < YmaxMV);
        check_property_default(gdist_subinterval_limit, size_t, parse_size, 10000)
        logger << "Creating MV gluon distribution with " << q2minMV << " < k2 < " << q2maxMV << ", " << YminMV << " < Y < " << YmaxMV << endl;
        gdist = new MVGluonDistribution(*satscale, lambdaMV, gammaMV, q2minMV, q2maxMV, YminMV, YmaxMV, gdist_subinterval_limit);
    }
    else if (gdist_type == "fMV") {
        double Ymin = min(Y);
        check_property_default(lambdaMV, double, parse_double, 0.241)
        check_property_default(gammaMV,  double, parse_double, 1)
        check_property_default(q2minMV,  double, parse_double, 1e-6)
        // q2max = (2 qxmax + sqrt(smax) / exp(Ymin))^2 + (2 qymax)^2
        check_property_default(q2maxMV,  double, parse_double, gsl_pow_2(2 * inf + sqs / exp(Ymin)) + gsl_pow_2(2 * inf))
        check_property(YMV, double, parse_double)
        logger << "Creating fMV gluon distribution with " << q2minMV << " < k2 < " << q2maxMV << ", Y = " << YMV << endl;
        assert(q2minMV < q2maxMV);
        check_property_default(gdist_subinterval_limit, size_t, parse_size, 10000)
        gdist = new FixedSaturationMVGluonDistribution(*satscale, lambdaMV, gammaMV, q2minMV, q2maxMV, YMV, gdist_subinterval_limit);
    }
    else if (gdist_type == "file") {
        check_property(gdist_position_filename, string, parse_string)
        check_property(gdist_momentum_filename, string, parse_string)
        check_property_default(xinit, double, parse_double, 0.01)
        logger << "Reading gluon distribution from " << gdist_position_filename << " (pos) and " << gdist_momentum_filename << " (mom)" << endl;
        gdist = new FileDataGluonDistribution(*satscale, gdist_position_filename, gdist_momentum_filename, xinit);
    }
    else {
        throw InvalidPropertyValueException<string>("gdist_type", gdist_type);
    }
    assert(gdist != NULL);
    if (trace_gdist) {
        logger << "Activating gluon distribution trace wrapper" << endl;
        gdist = new GluonDistributionTraceWrapper(gdist);
    }
    assert(gdist != NULL);

    // create coupling
    assert(cpl == NULL);
    check_property(coupling_type, string, parse_string)
    if (coupling_type == "fixed") {
        check_property(alphas, double, parse_double)
        cpl = new FixedCoupling(alphas);
    }
    else if (coupling_type == "running") {
        check_property(lambdaQCD, double, parse_double)
        check_property(regulator, double, parse_double)
        check_property_default(Ncbeta, double, parse_double, (11.0 * Nc - 2.0 * Nf) / 12.0)
        cpl = new LORunningCoupling(lambdaQCD, Ncbeta, regulator);
    }
    else {
        throw InvalidPropertyValueException<string>("coupling_type", coupling_type);
    }
    assert(cpl != NULL);
    
    contexts_created = true;
    // create contexts
    for (vector<double>::iterator pTit = pT.begin(); pTit != pT.end(); pTit++) {
        for (vector<double>::iterator Yit = Y.begin(); Yit != Y.end(); Yit++) {
            try {
                contexts.push_back(
                    Context(
                        x0,
                        A,
                        c,
                        lambda,
                        mu24pT2 ? 4*gsl_pow_2(*pTit) : mu2,
                        Nc,
                        Nf,
                        CF,
                        TR,
                        Sperp,
                        gsl_pow_2(*pTit),
                        sqs,
                        *Yit,
                        pdf_filename,
                        ff_filename,
                        projectile,
                        hadron,
                        integration_strategy,
                        cubature_iterations,
                        miser_iterations,
                        vegas_initial_iterations,
                        vegas_incremental_iterations,
                        quasi_iterations,
                        abserr,
                        relerr,
                        quasirandom_generator_type,
                        pseudorandom_generator_type,
                        pseudorandom_generator_seed,
                        gdist,
                        cpl));
            }
            catch (const InvalidKinematicsException& e) {
                logger << "Failed to create context at pT = " << *pTit << ", Y = " << *Yit << ": " << e.what() << endl;
            }
        }
    }
}

Context& ContextCollection::get_context(size_t n) {
    if (!contexts_created) {
        create_contexts();
    }
    return contexts[n];
}

Context& ContextCollection::operator[](size_t n) {
    if (!contexts_created) {
        create_contexts();
    }
    return get_context(n);
}

bool ContextCollection::empty() {
    return size() == 0;
}

size_t ContextCollection::size() {
    if (!contexts_created) {
        return options.count("pt") * options.count("y");
    }
    else {
        return contexts.size();
    }
}

void ContextCollection::set(string key, string value) {
    assert(!contexts_created); // TODO throw a proper exception here
    key = canonicalize(key);
    options.erase(key);
    options.insert(pair<string, string>(key, value));
}

void ContextCollection::erase(string key) {
    assert(!contexts_created);
    key = canonicalize(key);
    options.erase(key);
}

void ContextCollection::add(string key, string value) {
    assert(!contexts_created);
    key = canonicalize(key);
    // these are the keys that allow multiple values
    if (!(key == "pt" || key == "y")) {
        options.erase(key);
    }
    options.insert(pair<string, string>(key, value));
}

void ContextCollection::setup_defaults() {
    // This includes only those default values which are static
    options.insert(pair<string, string>("x0", "0.000304"));
    options.insert(pair<string, string>("lambda", "0.288"));
    options.insert(pair<string, string>("lambdaqcd", "0.24248711")); // sqrt(0.0588)
    options.insert(pair<string, string>("regulator", "1"));
    options.insert(pair<string, string>("mu2", "10"));
    options.insert(pair<string, string>("nc", "3"));
    options.insert(pair<string, string>("nf", "3"));
    options.insert(pair<string, string>("cf", "1.5"));
    options.insert(pair<string, string>("tr", "0.5"));
    options.insert(pair<string, string>("alphas", "0.2"));
    options.insert(pair<string, string>("pdf_filename", "mstw2008nlo.00.dat"));
    options.insert(pair<string, string>("ff_filename", "PINLO.DAT"));
    options.insert(pair<string, string>("integration_strategy", "VEGAS"));
    options.insert(pair<string, string>("cubature_iterations", "1000000"));
    options.insert(pair<string, string>("miser_iterations", "10000000"));
    options.insert(pair<string, string>("vegas_initial_iterations", "100000"));
    options.insert(pair<string, string>("vegas_incremental_iterations", "1000000"));
    options.insert(pair<string, string>("quasi_iterations", "1000000"));
    options.insert(pair<string, string>("abserr", "1e-20"));
    options.insert(pair<string, string>("relerr", "0"));
    // Adapted from the GSL source code - basically this reimplements gsl_rng_env_setup
    const char* qtype = getenv("GSL_QRNG_TYPE");
    options.insert(pair<string, string>("quasirandom_generator_type", qtype == NULL ? "halton" : qtype));
    const char* type = getenv("GSL_RNG_TYPE");
    options.insert(pair<string, string>("pseudorandom_generator_type", type == NULL ? "mt19937" : type));
    const char* seed = getenv("GSL_RNG_SEED");
    options.insert(pair<string, string>("pseudorandom_generator_seed", seed == NULL ? "0" : seed));
}

void ContextCollection::read_config(istream& in) {
    string line;
    do {
        getline(in, line);
        if (line.size() > 2 && line[0] != '#') {
            // Split the line into two pieces on the '=' character
            // The first piece becomes the key, the second becomes the value
            vector<string> kv = split(line, "\n=", 2);
            string key = canonicalize(kv[0]);
            // split the value on commas
            vector<string> v = split(kv[1], ",");
            for (vector<string>::iterator it = v.begin(); it != v.end(); it++) {
                string value = trim(*it, " \n\t");
                add(key, value);
            }
        }
    } while (!in.eof());
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
}

std::ostream& operator<<(std::ostream& out, Context& ctx) {
    out << "x0\t= "         << ctx.x0           << endl;
    out << "A\t= "          << ctx.A            << endl;
    out << "c\t= "          << ctx.c            << endl;
    out << "lambda\t= "     << ctx.lambda       << endl;
    out << "mu2\t= "        << ctx.mu2          << endl;
    out << "Nc\t= "         << ctx.Nc           << endl;
    out << "Nf\t= "         << ctx.Nf           << endl;
    out << "CF\t= "         << ctx.CF           << endl;
    out << "TR\t= "         << ctx.TR           << endl;
    out << "Sperp\t= "      << ctx.Sperp        << endl;
    out << "pT2\t= "        << ctx.pT2          << endl;
    out << "sqs\t= "        << ctx.sqs          << endl;
    out << "Y\t= "          << ctx.Y            << endl;
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
    out << "gluon distribution\t = " << ctx.gdist << endl;
    out << "coupling\t = " << ctx.cpl << endl;
    out << "quasirandom generator type: " <<  ctx.quasirandom_generator_type->name << endl;
    out << "pseudorandom generator type: " <<  ctx.pseudorandom_generator_type->name << endl;
    out << "pseudorandom generator seed: " << ctx.pseudorandom_generator_seed << endl;
    return out;
}

template<typename T>
ostream& operator<<(ostream& out, vector<T>& vec) {
    for (typename vector<T>::iterator it = vec.begin(); it != vec.end(); it++) {
        if (it != vec.begin()) {
            out << ", ";
        }
        out << *it;
    }
    return out;
}

std::istream& operator>>(std::istream& in, ContextCollection& cc) {
    cc.read_config(in);
    return in;
}

std::ostream& operator<<(std::ostream& out, ContextCollection& cc) {
    string last_key;
    for (multimap<string, string>::iterator it = cc.options.begin(); it != cc.options.end(); it++) {
        if (last_key == it->first) {
            out << ", " << it->second;
        }
        else {
            if (!last_key.empty()) {
                out << endl;
            }
            out << it->first << " = " << it->second;
        }
        last_key = it->first;
    }
    out << endl;
    return out;
}

ContextCollection::iterator ContextCollection::begin() {
    if (contexts.empty()) {
        create_contexts();
    }
    return contexts.begin();
}
ContextCollection::iterator ContextCollection::end() {
    if (contexts.empty()) {
        create_contexts();
    }
    return contexts.end();
}


#ifdef CONTEXT_TEST
const double inf = 10;
ostream& logger = cerr;

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <filename.cfg>" << endl;
        return 1;
    }
    ContextCollection cc(argv[1]);
    try {
        cc.create_contexts();
    }
    catch (const exception& e) {
        cout << "Error in parsing: " << e.what() << endl;
        return 1;
    }
    cout << "Successfully parsed " << argv[1];
    if (cc.empty()) {
        cout << endl << cc << "No contexts defined!" << endl;
    }
    else {
        cout << " into " << cc.size() << " contexts" << endl;
        cout << cc;
        for (ContextCollection::iterator it = cc.begin(); it != cc.end(); it++) {
            cout << "-----------------------------------" << endl << "Context:" << endl << *it;
        }
    }
    return 0;
}
#endif

