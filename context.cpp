#include <fstream>
#include <iostream>
#include <map>
#include <cstdlib>
#include "context.h"
#include "utils.h"

using namespace std;

Context ContextCollection::get_context(size_t n) {
#define check_property(p) if (p == unset) { cerr << "No value for " #p << endl; err = true; }
#define check_container(p) if (p.empty()) { cerr << "No value for " #p << endl; err = true; }
    bool err = false;
    check_property(x0)
    check_property(A)
    check_property(c)
    check_property(lambda)
    check_property(lambdaQCD)
    check_property(mu2)
    check_property(Nc)
    check_property(Nf)
    check_property(CF)
    check_property(TR)
    check_property(Sperp)
    check_container(pT2)
    check_property(sqs)
    check_container(Y)
    check_container(pdf_filename)
    check_container(ff_filename)
    if (err) {
        exit(1);
    }
#undef check_property
#undef check_container
    double l_pT2 = pT2[n / Y.size()];
    double l_Y = Y[n % Y.size()];
    Context ctx(
        x0,
        A,
        c,
        lambda,
        lambdaQCD,
        mu2,
        Nc,
        Nf,
        CF,
        TR,
        Sperp,
        l_pT2,
        sqs,
        l_Y,
        pdf_filename,
        ff_filename,
        gdist,
        cpl);
    return ctx;
}

const string canonicalize(const string& key) {
    if (key == "lambda_QCD" || key == "LAMBDA_QCD") {
        return "lambdaQCD";
    }
    else if (key == "mu^2") {
        return "mu2";
    }
    else if (key == "ncolors" || key == "nc") {
        return "NC";
    }
    else if (key == "nflavors" || key == "nf") {
        return "NF";
    }
    else if (key == "colorfactor" || key == "cf") {
        return "CF";
    }
    else if (key == "tr") {
        return "TR";
    }
    else if (key == "alphas_bar" || key == "alpha_s_bar") {
        return "alphasbar";
    }
    else {
        return key;
    }
}


void ContextCollection::read_config(istream& in) {
    multimap<string, string> options;

    // read options
    string line;
    do {
        getline(in, line);
        if (line.size() > 2 && line[0] != '#') {
            // Split the line into two pieces on the '=' character
            // The first piece becomes the key, the second becomes the value
            vector<string> kv = split(line, " \n\t=", 2);
            const string key = canonicalize(kv[0]);
            // split the value on commas
            vector<string> v = split(kv[1], ",");
            for (vector<string>::iterator it = v.begin(); it != v.end(); it++) {
                pair<string, string> p(key, trim(*it, " \n\t"));
                options.insert(p);
            }
        }
    } while (!in.eof());

    // now process the options
    for (map<string, string>::const_iterator it = options.begin(); it != options.end(); it++) {
        const string key = it->first;
        const string value = it->second;
        if (key == "x0") {
            x0 = atof(value.c_str());
        }
        else if (key == "A") {
            A = atoi(value.c_str());
        }
        else if (key == "c") {
            c = atof(value.c_str());
        }
        else if (key == "lambda") {
            lambda = atof(value.c_str());
        }
        else if (key == "lambdaQCD") {
            lambdaQCD = atof(value.c_str());
        }
        else if (key == "mu2") {
            mu2 = atof(value.c_str());
        }
        else if (key == "NC") {
            Nc = atoi(value.c_str());
        }
        else if (key == "NF") {
            Nf = atoi(value.c_str());
        }
        else if (key == "CF") {
            CF = atof(value.c_str());
        }
        else if (key == "TR") {
            TR = atof(value.c_str());
        }
        else if (key == "Sperp") {
            Sperp = atof(value.c_str());
        }
        else if (key == "pT2") {
            pT2.push_back(atof(value.c_str()));
        }
        else if (key == "sqs") {
            sqs = atof(value.c_str());
        }
        else if (key == "Y") {
            Y.push_back(atof(value.c_str()));
        }
        else if (key == "pdf_filename") {
            pdf_filename = value;
        }
        else if (key == "ff_filename") {
            ff_filename = value;
        }
        else if (key == "coupling") {
            if (value == "fixed") {
                multimap<string, string>::const_iterator subit = options.find("alphasbar");
                if (subit == options.end()) {
                    cerr << "No value for fixed coupling" << endl;
                    continue;
                }
                double alphasbar = atof(subit->second.c_str());
                cpl = new FixedCoupling(alphasbar);
            }
            else if (value == "running") {
                multimap<string, string>::const_iterator subit = options.find("lambdaQCD");
                if (subit == options.end()) {
                    cerr << "No value for lambdaQCD";
                    continue;
                }
                double lambdaQCD = atof(subit->second.c_str());
                subit = options.find("beta");
                if (subit == options.end()) {
                    cerr << "No value for beta";
                    continue;
                }
                double beta = atof(subit->second.c_str());
                subit = options.find("regulator");
                if (subit == options.end()) {
                    cerr << "No value for regulator";
                    continue;
                }
                double regulator = atof(subit->second.c_str());
                cpl = new LORunningCoupling(lambdaQCD, beta, regulator);
            }
        }
        else if (key == "gluon distribution" || key == "gdist") {
            if (value == "GBW") {
                gdist = new GBWGluonDistribution();
            }
            else if (value == "MV") {
                multimap<string, string>::const_iterator subit = options.find("lambdaMV");
                if (subit == options.end()) {
                    cerr << "No value for Lambda MV" << endl;
                    continue;
                }
                double lambdaMV = atof(subit->second.c_str());
                subit = options.find("q2min");
                if (subit == options.end()) {
                    cerr << "No value for q2min" << endl;
                    continue;
                }
                double q2min = atof(subit->second.c_str());
                subit = options.find("q2max");
                if (subit == options.end()) {
                    cerr << "No value for q2max" << endl;
                    continue;
                }
                double q2max = atof(subit->second.c_str());
                subit = options.find("Qs2min");
                if (subit == options.end()) {
                    cerr << "No value for Qs2min" << endl;
                    continue;
                }
                double Qs2min = atof(subit->second.c_str());
                subit = options.find("Qs2max");
                if (subit == options.end()) {
                    cerr << "No value for Qs2max" << endl;
                    continue;
                }
                double Qs2max = atof(subit->second.c_str());
                gdist = new MVGluonDistribution(lambdaMV, q2min, q2max, Qs2min, Qs2max);
            }
        }
        else {
            cerr << "unrecognized property " << key << endl;
        }
    }
}

std::ostream& operator<<(std::ostream& out, Context& ctx) {
    out << "x0\t= "         << ctx.x0           << endl;
    out << "A\t= "          << ctx.A            << endl;
    out << "c\t= "          << ctx.c            << endl;
    out << "lambda\t= "     << ctx.lambda       << endl;
    out << "lambdaQCD\t= "  << ctx.lambdaQCD    << endl;
    out << "mu2\t= "        << ctx.mu2          << endl;
    out << "Nc\t= "         << ctx.Nc           << endl;
    out << "Nf\t= "         << ctx.Nf           << endl;
    out << "CF\t= "         << ctx.CF           << endl;
    out << "TR\t= "         << ctx.TR           << endl;
    out << "Sperp\t= "      << ctx.Sperp        << endl;
    out << "pT2\t= "        << ctx.pT2          << endl;
    out << "sqs\t= "        << ctx.sqs          << endl;
    out << "Y\t= "          << ctx.Y            << endl;
    out <<"pdf_filename\t= "<< ctx.pdf_filename << endl;
    out << "ff_filename\t= "<< ctx.ff_filename  << endl;
//     out << "gdist\t="       << ctx.gdist        << endl;
//     out << "coupling\t="    << ctx.cpl          << endl;
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
    out << "x0           = " << cc.x0           << endl;
    out << "A            = " << cc.A            << endl;
    out << "c            = " << cc.c            << endl;
    out << "lambda       = " << cc.lambda       << endl;
    out << "lambdaQCD    = " << cc.lambdaQCD    << endl;
    out << "mu2          = " << cc.mu2          << endl;
    out << "Nc           = " << cc.Nc           << endl;
    out << "Nf           = " << cc.Nf           << endl;
    out << "CF           = " << cc.CF           << endl;
    out << "TR           = " << cc.TR           << endl;
    out << "Sperp        = " << cc.Sperp        << endl;
    out << "pT2          = " << cc.pT2          << endl;
    out << "sqs          = " << cc.sqs          << endl;
    out << "Y            = " << cc.Y            << endl;
    out << "pdf_filename = " << cc.pdf_filename << endl;
    out << "ff_filename  = " << cc.ff_filename  << endl;
//     out << "gdist\t="       << cc.gdist        << endl;
//     out << "coupling\t="    << cc.cpl          << endl;
    return out;
}

ContextCollection::iterator ContextCollection::begin() {
    return ContextCollectionIterator(*this, 0);
}
ContextCollection::iterator ContextCollection::end() {
    return ContextCollectionIterator(*this, size());
}


#ifdef CONTEXT_TEST
int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <filename.cfg>" << endl;
        return 1;
    }
    ContextCollection cc(argv[1]);
    cout << "Successfully parsed " << argv[1] << " into " << cc.size() << " contexts" << endl;
    cout << cc;
}
#endif

