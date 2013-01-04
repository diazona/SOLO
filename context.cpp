#include <cassert>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <gsl/gsl_math.h>
#include "context.h"
#include "utils.h"

#define check_property_default(p, typename, parse, default) \
    itit = options.equal_range(canonicalize(#p));\
    typename p;\
    if (itit.first == itit.second) {\
        p = default;\
    }\
    else {\
        p = parse(itit);\
    }
#define check_property(p, typename, parse) \
    itit = options.equal_range(canonicalize(#p));\
    typename p;\
    if (itit.first == itit.second) {\
        cerr << "No value for " #p << endl;\
        err = true;\
    }\
    else {\
        p = parse(itit);\
    }

using namespace std;

extern const double inf;

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
    else if (key == "alphas_bar" || key == "alpha_s_bar") {
        return "alphasbar";
    }
    else if (key == "coupling" || key == "cpl") {
        return "coupling_type";
    }
    else if (key == "gdist" || key == "gluon_distribution" || key == "gluon distribution" || key == "gluon dist") {
        return "gdist_type";
    }
    else {
        return key;
    }
}

size_t parse_size(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    multimap<string, string>::iterator el = range.first;
    assert(++el == range.second);
    return (size_t)strtoul(range.first->second.c_str(), NULL, 0);
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

vector<double> parse_vector(pair<multimap<string, string>::iterator, multimap<string, string>::iterator> range) {
    vector<double> v;
    for (multimap<string, string>::iterator it = range.first; it != range.second; it++) {
        v.push_back(strtod(it->second.c_str(), NULL));
    }
    return v;
}

void ContextCollection::create_contexts() {
    pair<multimap<string, string>::iterator, multimap<string, string>::iterator> itit;
    bool err = false;

    check_property(x0,           double, parse_double)
    check_property(A,            double, parse_double)
    check_property(c,            double, parse_double)
    check_property(lambda,       double, parse_double)
    check_property(mu2,          double, parse_double)
    check_property(Nc,           double, parse_double)
    check_property(Nf,           double, parse_double)
    check_property(CF,           double, parse_double)
    check_property(TR,           double, parse_double)
    check_property(Sperp,        double, parse_double)
    check_property(sqs,          double, parse_double)
    check_property(pT,           vector<double>, parse_vector)
    check_property(Y,            vector<double>, parse_vector)
    check_property(pdf_filename, string, parse_string)
    check_property(ff_filename,  string, parse_string)
    check_property(miser_iterations, size_t, parse_size)
    check_property(vegas_initial_iterations, size_t, parse_size)
    check_property(vegas_incremental_iterations, size_t, parse_size)
    if (err) {
        exit(1);
    }
    
    // create gluon distribution
    if (gdist == NULL) {
        check_property(gdist_type, string, parse_string)
        if (err) {
            exit(1);
        }
        if (gdist_type == "GBW") {
            gdist = new GBWGluonDistribution();
        }
        else if (gdist_type == "MV") {
            double Ymin = min(Y);
            double Ymax = max(Y);
            double pTmin = min(pT);
            double pTmax = max(pT);
            double Q02x0lambda = c * pow(A, 1./3.) * pow(x0, lambda);
            check_property(lambdaMV, double, parse_double)
            check_property_default(q2min,  double, parse_double, 1e-6)
            // q2max = (2 qxmax + sqrt(smax) / exp(Ymin))^2 + (2 qymax)^2
            check_property_default(q2max,  double, parse_double, gsl_pow_2(2 * inf + sqs / exp(Ymin)) + gsl_pow_2(2 * inf))
            // Qs2min = c A^1/3 Q02 (x0 / exp(-2Ymin))^λ
            check_property_default(Qs2min, double, parse_double, Q02x0lambda * exp(2 * lambda * Ymin))
            // Qs2max = c A^1/3 Q02 x0^λ / (pT / sqs * exp(-Ymin))^λ
            check_property_default(Qs2max, double, parse_double, Q02x0lambda * pow(pTmin / sqs * exp(-Ymax), -lambda))
            cerr << "Creating MV gluon distribution with " << q2min << " < k2 < " << q2max << ", " << Qs2min << " < Qs2 < " << Qs2max << endl;
            assert(q2min < q2max);
            assert(Qs2min < Qs2max);
            gdist = new MVGluonDistribution(lambdaMV, q2min, q2max, Qs2min, Qs2max);
        }
        else {
            cerr << "Invalid value '" << gdist_type << "' for gdist_type!" << endl;
            err = true;
        }
        if (err) {
            exit(1);
        }
    }

    // create coupling
    if (cpl == NULL) {
        check_property(coupling_type, string, parse_string)
        if (err) {
            exit(1);
        }
        if (coupling_type == "fixed") {
            check_property(alphasbar, double, parse_double)
            cpl = new FixedCoupling(alphasbar);
        }
        else if (coupling_type == "running") {
            check_property(lambdaQCD, double, parse_double)
            check_property(beta,      double, parse_double)
            check_property(regulator, double, parse_double)
            cpl = new LORunningCoupling(lambdaQCD, beta, regulator);
        }
        else {
            cerr << "Invalid value '" << coupling_type << "' for coupling_type!" << endl;
            err = true;
        }
        if (err) {
            exit(1);
        }
    }
    
    // create contexts
    for (vector<double>::iterator pTit = pT.begin(); pTit != pT.end(); pTit++) {
        for (vector<double>::iterator Yit = Y.begin(); Yit != Y.end(); Yit++) {
            contexts.push_back(Context(x0, A, c, lambda, mu2, Nc, Nf, CF, TR, Sperp, gsl_pow_2(*pTit), sqs, *Yit, pdf_filename, ff_filename, miser_iterations, vegas_initial_iterations, vegas_incremental_iterations, gdist, cpl));
        }
    }
}

Context& ContextCollection::get_context(size_t n) {
    if (contexts.empty()) {
        create_contexts();
    }
    return contexts[n];
}

Context& ContextCollection::operator[](size_t n) {
    if (contexts.empty()) {
        create_contexts();
    }
    return get_context(n);
}

bool ContextCollection::empty() {
    return size() == 0;
}

size_t ContextCollection::size() {
    if (contexts.empty()) {
        return options.count("pt") * options.count("y");
    }
    else {
        return contexts.size();
    }
}

void ContextCollection::set(string key, string value) {
    assert(contexts.empty());
    key = canonicalize(key);
    options.erase(key);
    options.insert(pair<string, string>(key, value));
}

void ContextCollection::erase(string key) {
    assert(contexts.empty());
    key = canonicalize(key);
    options.erase(key);
}

void ContextCollection::add(string key, string value) {
    assert(contexts.empty());
    key = canonicalize(key);
    // these are the keys that allow multiple values
    if (!(key == "pt" || key == "y")) {
        options.erase(key);
    }
    options.insert(pair<string, string>(key, value));
}

void ContextCollection::setup_defaults() {
    options.insert(pair<string, string>("x0", "0.000304"));
    options.insert(pair<string, string>("lambda", "0.288"));
    options.insert(pair<string, string>("lambdamv", "0.24"));
    options.insert(pair<string, string>("lambdaqcd", "0.24248711")); // sqrt(0.0588)
    options.insert(pair<string, string>("mu2", "10"));
    options.insert(pair<string, string>("nc", "3"));
    options.insert(pair<string, string>("nf", "3"));
    options.insert(pair<string, string>("cf", "1.5"));
    options.insert(pair<string, string>("tr", "0.5"));
    options.insert(pair<string, string>("pdf_filename", "mstw2008nlo.00.dat"));
    options.insert(pair<string, string>("ff_filename", "PINLO.DAT"));
    options.insert(pair<string, string>("miser_iterations", "10000000"));
    options.insert(pair<string, string>("vegas_initial_iterations", "100000"));
    options.insert(pair<string, string>("vegas_incremental_iterations", "1000000"));
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
    out << "miser_iterations\t= " << ctx.miser_iterations << endl;
    out << "vegas_initial_iterations\t= " << ctx.vegas_initial_iterations << endl;
    out << "vegas_incremental_iterations\t= " << ctx.vegas_incremental_iterations << endl;
    out << "gluon distribution\t = " << ctx.gdist << endl;
    out << "coupling\t = " << ctx.cpl << endl;
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

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <filename.cfg>" << endl;
        return 1;
    }
    ContextCollection cc(argv[1]);
    cout << "Successfully parsed " << argv[1];
    if (cc.empty()) {
        cout << endl << cc << "No contexts defined!" << endl;
    }
    else {
        cout << " into " << cc.size() << " contexts" << endl;
        cout << cc;
    }
}
#endif

