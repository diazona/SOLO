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

using namespace std;

const int SUCCESS = 0;

static bool trace = false;
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

class ResultsCalculator {
private:
    ContextCollection& cc;
    ThreadLocalContext& tlctx;
    vector<HardFactorList*> hfgroups;
    double* real;
    double* imag;
    double* error;
public:
    ResultsCalculator(ContextCollection& cc, ThreadLocalContext& tlctx, vector<HardFactorList*> hfgroups) : cc(cc), tlctx(tlctx), hfgroups(hfgroups) {
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
    void result(size_t ccindex, size_t hfindex, double* real, double* imag, double* error);
    void calculate();
};

void ResultsCalculator::result(size_t ccindex, size_t hfindex, double* real, double* imag, double* error) {
    size_t index = index_from(ccindex, hfindex);
    (*real) = this->real[index];
    (*imag) = this->imag[index];
    (*error) = this->error[index];
}

void ResultsCalculator::calculate() {
    double* l_real = real;
    double* l_imag = imag;
    double* l_error = error;
    for (ContextCollection::iterator it = cc.begin(); it != cc.end(); it++) {
        Context ctx = *it;
        cerr << "Beginning calculation at pT = " << sqrt(ctx.pT2) << ", Y = " << ctx.Y << endl;
        for (vector<HardFactorList*>::iterator hit = hfgroups.begin(); hit != hfgroups.end(); hit++) {
            Integrator* integrator = new Integrator(&ctx, &tlctx, strategy, **hit);
            if (trace) {
                integrator->set_callback(write_data_point);
            }
            else if (minmax) {
                integrator->set_callback(store_minmax);
            }
            integrator->integrate(l_real++, l_imag++, l_error++);
            delete integrator;
        }
        cerr << "...done" << endl;
    }
}

HardFactorList* parse_hf_spec(const string& spec) {
    vector<string> hfnames = split(spec, ", ");
    HardFactorList* hfobjs = new HardFactorList();
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
            cerr << "No such hard factor " << *it << endl;
            exit(1);
        }
        hfobjs->push_back(hf);
    }
    assert(hfobjs != NULL);
    assert(hfobjs->size() > 0);
    return hfobjs;
}

static const char* default_lo_spec = "p:h02qq,p:h02gg";
static const char* default_nlo_spec = "p:h12qq,p:h14qq,p:h12gg,m:h12qqbar,m:h16gg,p:h112gq,p:h122gq,m:h14gq,p:h112qg,p:h122qg,m:h14qg";

int main(int argc, char** argv) {
    bool separate = false;
    enum {GBW, MV} gdist_type = GBW; 
    Coupling* cpl = new FixedCoupling(0.2 / (2*M_PI));
    vector<HardFactorList*> hfgroups;
    vector<string> hfgnames;
    vector<double> pT2;
    vector<double> Y;
    ContextCollection cc;
    ifstream config;
    bool config_is_read = false;

    // process options
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--miser")==0) {
            strategy = MC_MISER;
        }
        else if (strcmp(argv[i], "--vegas")==0) {
            strategy = MC_VEGAS;
        }
        else if (strcmp(argv[i], "--trace")==0) {
            trace = true;
        }
        else if (strcmp(argv[i], "--minmax")==0) {
            minmax = true;
        }
        else if (strcmp(argv[i], "--separate")==0) {
            separate = true;
        }
        else if (strcmp(argv[i], "MV") == 0) {
            gdist_type = MV;
        }
        else if (strcmp(argv[i], "GBW") == 0) {
            gdist_type = GBW;
        }
        else if (argv[i][0] == 'h' || argv[i][1] == ':') {
            hfgroups.push_back(parse_hf_spec(argv[i]));
            hfgnames.push_back(argv[i]);
        }
        else if (strcmp(argv[i], "lo") == 0) {
            hfgroups.push_back(parse_hf_spec(default_lo_spec));
            hfgnames.push_back(argv[i]);
        }
        else if (strcmp(argv[i], "nlo") == 0) {
            hfgroups.push_back(parse_hf_spec(default_nlo_spec));
            hfgnames.push_back(argv[i]);
        }
        else if (::isdigit(argv[i][0])) {
            vector<string> pTnums = split(argv[i], ", ");
            for (vector<string>::iterator it = pTnums.begin(); it != pTnums.end(); it++) {
                double d = strtod(it->c_str(), NULL);
                if (errno != 0) {
                    cerr << "Error " << errno << " parsing " << argv[i] << " as floating-point" << endl;
                    cerr << ERANGE << ": range error" << endl;
                    exit(1);
                }
                pT2.push_back(d*d);
            }
        }
        else {
            // try opening as a file
            config.open(argv[i]);
            if (config.good()) {
                if (config_is_read) {
                    cerr << "Warning: reading extra config file " << argv[i] << endl;
                }
                config >> cc;
                config.close();
                config_is_read = true;
            }
            else {
                cerr << "Unrecognized argument " << argv[i] << endl;
                exit(1);
            }
        }
    }

    if (!config_is_read) {
        cerr << "No config file specified! (one argument must name a config file)" << endl;
        exit(1);
    }
    if (pT2.empty() && cc.pT2.empty()) {
        cerr << "No momenta specified!" << endl;
        exit(1);
    }
    if (cc.Y.empty()) {
        cerr << "No rapidities specified!" << endl;
        exit(1);
    }
    if (!pT2.empty()) {
        cc.pT2 = pT2;
    }
    if (hfgroups.empty()) {
        hfgroups.push_back(parse_hf_spec(default_lo_spec));
        hfgnames.push_back("lo");
        hfgroups.push_back(parse_hf_spec(default_nlo_spec));
        hfgnames.push_back("nlo");
    }
    assert(!hfgroups.empty());
    assert(hfgroups.size() == hfgnames.size());
    
    gsl_rng_env_setup();
    cc.A = 197;
    cc.c = 0.56;
    cc.mu2 = 10;
    cc.Sperp = 1.0;
    cc.sqs = 200;
    cc.cpl = cpl;
    ThreadLocalContext tlctx(cc);
    double pT2min = min(cc.pT2);
    double pT2max = max(cc.pT2);
    double Ymin = min(cc.Y);
    double Ymax = max(cc.Y);
    double k2min, k2max, Qs2min, Qs2max;
    switch (gdist_type) {
        case MV:
            // TODO: check usage of Ymax and Ymin in these formulas
            k2min = 1e-6;
            k2max = gsl_pow_2(2 * inf + cc.sqs / exp(Ymin)) + gsl_pow_2(2 * inf); // (2 qxmax + sqrt(smax) / exp(Ymin))^2 + (2 qymax)^2
            Qs2min = cc.Q02x0lambda() * exp(2 * cc.lambda * Ymin); // c A^1/3 Q02 (x0 / exp(-2Ymin))^λ
            Qs2max = cc.Q02x0lambda() * pow(sqrt(pT2max) / cc.sqs * exp(-Ymin), -cc.lambda); // c A^1/3 Q02 x0^λ / (pT / sqs * exp(-Ymin))^λ
            cerr << "Creating MV gluon distribution with " << k2min << " < k2 < " << k2max << ", " << Qs2min << " < Qs2 < " << Qs2max << endl;
            assert(k2min < k2max);
            assert(Qs2min < Qs2max);
            cc.gdist = new MVGluonDistribution(
                0.24,  // TODO: replace with LambdaMV parameter
                k2min, // k2min
                k2max, // k2max
                Qs2min,// Qs2min
                Qs2max // Qs2max
            );
            break;
        case GBW:
        default:
            cc.gdist = new GBWGluonDistribution();
            break;
    }

    ResultsCalculator* rc = new ResultsCalculator(cc, tlctx, hfgroups);
    rc->calculate();
    cout << "pT\tY\t";
    for (vector<string>::iterator it = hfgnames.begin(); it != hfgnames.end(); it++) {
        cout << *it << "\t";
    }
    cout << "total" << endl;
    double l_real, l_imag, l_error;
    for (size_t ccindex = 0; ccindex < cc.size(); ccindex++) {
        cout << sqrt(cc[ccindex].pT2) << "\t";
        cout << sqrt(cc[ccindex].Y) << "\t";

        double total = 0;
        for (size_t hfgindex = 0; hfgindex < hfgroups.size(); hfgindex++) {
            rc->result(ccindex, hfgindex, &l_real, &l_imag, &l_error);
            cout << l_real << "±" << l_error << "\t";
            total += l_real;
        }
        cout << total << endl;
    }
    
    if (minmax) {
        cerr << "xx\t" << min_ictx.xx << "\t" << max_ictx.xx << "\t" << endl;
        cerr << "xy\t" << min_ictx.xy << "\t" << max_ictx.xy << "\t" << endl;
        cerr << "yx\t" << min_ictx.yx << "\t" << max_ictx.yx << "\t" << endl;
        cerr << "yy\t" << min_ictx.yy << "\t" << max_ictx.yy << "\t" << endl;
        cerr << "bx\t" << min_ictx.bx << "\t" << max_ictx.bx << "\t" << endl;
        cerr << "by\t" << min_ictx.by << "\t" << max_ictx.by << "\t" << endl;
        cerr << "q1x\t" << min_ictx.q1x << "\t" << max_ictx.q1x << "\t" << endl;
        cerr << "q1y\t" << min_ictx.q1y << "\t" << max_ictx.q1y << "\t" << endl;
        cerr << "q2x\t" << min_ictx.q2x << "\t" << max_ictx.q2x << "\t" << endl;
        cerr << "q2y\t" << min_ictx.q2y << "\t" << max_ictx.q2y << "\t" << endl;
        cerr << "q3x\t" << min_ictx.q3x << "\t" << max_ictx.q3x << "\t" << endl;
        cerr << "q3y\t" << min_ictx.q3y << "\t" << max_ictx.q3y << "\t" << endl;
        cerr << "z\t" << min_ictx.z << "\t" << max_ictx.z << "\t" << endl;
        cerr << "xi\t" << min_ictx.xi << "\t" << max_ictx.xi << "\t" << endl;
        cerr << "xip\t" << min_ictx.xiprime << "\t" << max_ictx.xiprime << "\t" << endl;
        cerr << "kT\t" << min_ictx.kT << "\t" << max_ictx.kT << "\t" << endl;
        cerr << "kT2\t" << min_ictx.kT2 << "\t" << max_ictx.kT2 << "\t" << endl;
        cerr << "xp\t" << min_ictx.xp << "\t" << max_ictx.xp << "\t" << endl;
        cerr << "xg\t" << min_ictx.xg << "\t" << max_ictx.xg << "\t" << endl;
        cerr << "Qs2\t" << min_ictx.Qs2 << "\t" << max_ictx.Qs2 << "\t" << endl;
        cerr << "alphasbar\t" << min_ictx.alphasbar << "\t" << max_ictx.alphasbar << "\t" << endl;
    }
    
    for (vector<HardFactorList*>::iterator it = hfgroups.begin(); it != hfgroups.end(); it++) {
        delete *it;
    }
    delete cc.gdist;
    delete cpl;
    return 0;
}
