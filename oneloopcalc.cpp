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

static IntegrationContext min_ictx(NULL);
static IntegrationContext max_ictx(NULL);

#define store(property) \
  min_ictx.property = min_ictx.property == 0 ? ictx->property : min(min_ictx.property, ictx->property); \
  max_ictx.property = max_ictx.property == 0 ? ictx->property : max(max_ictx.property, ictx->property);

void store_minmax(const IntegrationContext* ictx, const double real, const double imag) {
    if (ictx == NULL) {
        return;
    }
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
    store(xi);
    store(z);
    store(xiprime);
}

void write_nonzero(const IntegrationContext* ictx, const double real, const double imag) {
    if (real != 0 || imag != 0) {
        cerr << real << "\t" << imag << endl;
    }
}

class ResultsCalculator {
private:
    Context* ctx;
    vector<double> pTlist;
    vector<HardFactorList*> hfgroups;
    double* real;
    double* imag;
    double* error;
public:
    ResultsCalculator(Context* ctx, vector<double> pTlist, vector<HardFactorList*> hfgroups) : ctx(ctx), pTlist(pTlist), hfgroups(hfgroups) {
        size_t hflen = hfgroups.size();
        size_t pTlen = pTlist.size();
        real = new double[hflen * pTlen];
        imag = new double[hflen * pTlen];
        error = new double[hflen * pTlen];
    };
    ~ResultsCalculator() {
        delete[] real;
        delete[] imag;
        delete[] error;
    };
    size_t index_from(size_t pTindex, size_t hfindex) {
        return pTindex * hfgroups.size() + hfindex;
    };
    void result(size_t pTindex, size_t hfindex, double* real, double* imag, double* error);
    void calculate();
};

void ResultsCalculator::result(size_t pTindex, size_t hfindex, double* real, double* imag, double* error) {
    size_t index = index_from(pTindex, hfindex);
    (*real) = this->real[index];
    (*imag) = this->imag[index];
    (*error) = this->error[index];
}

void ResultsCalculator::calculate() {
    double* l_real = real;
    double* l_imag = imag;
    double* l_error = error;
    for (vector<double>::iterator it = pTlist.begin(); it != pTlist.end(); it++) {
        double pT = *it;
        ctx->pT2 = pT*pT;
        ctx->recalculate();
        cerr << "Beginning calculation at pT = " << pT << endl;
        for (vector<HardFactorList*>::iterator hit = hfgroups.begin(); hit != hfgroups.end(); hit++) {
            Integrator* integrator = new Integrator(ctx, strategy, **hit);
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

// from http://stackoverflow.com/a/236803/56541
vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while(getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

HardFactorList* parse_hf_spec(const string& spec) {
    vector<string> hfnames;
    split(spec, ',', hfnames);
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
    vector<double> pT;

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
        else if (argv[i][0] == 'h') {
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
            vector<string> pTnums;
            split(argv[i], ',', pTnums);
            for (vector<string>::iterator it = pTnums.begin(); it != pTnums.end(); it++) {
                double d = strtod(it->c_str(), NULL);
                if (errno != 0) {
                    cerr << "Error " << errno << " parsing " << argv[i] << " as floating-point" << endl;
                    cerr << ERANGE << ": range error" << endl;
                    exit(1);
                }
                pT.push_back(d);
            }
        }
        else {
            cerr << "Unrecognized argument " << argv[i] << endl;
        }
    }

    if (pT.size() == 0) {
        cerr << "No momenta specified!" << endl;
        exit(1);
    }
    if (hfgroups.size() == 0) {
        hfgroups.push_back(parse_hf_spec(default_lo_spec));
        hfgnames.push_back("lo");
        hfgroups.push_back(parse_hf_spec(default_nlo_spec));
        hfgnames.push_back("nlo");
    }
    assert(hfgroups.size() > 0);
    assert(hfgroups.size() == hfgnames.size());
    
    gsl_rng_env_setup();
    Context gctx(
      0.000304, // x0
      197,      // A
      0.56,     // c
      0.288,    // lambda
      10,       // mu2
      3,        // Nc
      3,        // Nf
      1.5,      // CF
      0.5,      // TR
      1.0,      // Sperp
      1.0,      // pT2 (dummy value)
      200,      // sqs
      3.2,      // Y
      NULL,     // to be inserted later
      cpl,
      "mstw2008nlo.00.dat", "PINLO.DAT");
    double k2min, k2max, Qs2min, Qs2max;
    switch (gdist_type) {
        case MV:
            k2min = 1e-12;
            k2max = gsl_pow_2(inf + gctx.sqs / exp(gctx.Y)); // qmax + sqrt(smax) / exp(Ymin)
            Qs2min = gctx.Q02x0lambda; // c A^1/3 Q02 x0^λ
            Qs2max = gctx.Q02x0lambda / exp(-2 * gctx.lambda * gctx.Y); // c A^1/3 Q02 (x0 / exp(-2Y))^λ
            assert(k2min < k2max);
            assert(Qs2min < Qs2max);
            gctx.gdist = new MVGluonDistribution(
                0.24,  // TODO: replace with LambdaMV parameter
                k2min, // k2min
                k2max, // k2max
                Qs2min,// Qs2min
                Qs2max // Qs2max
            );
            break;
        case GBW:
        default:
            gctx.gdist = new GBWGluonDistribution();
            break;
    }

    ResultsCalculator* rc = new ResultsCalculator(&gctx, pT, hfgroups);
    rc->calculate();
    cout << "pT\t";
    for (vector<string>::iterator it = hfgnames.begin(); it != hfgnames.end(); it++) {
        cout << *it << "\t";
    }
    cout << "total" << endl;
    double l_real, l_imag, l_error;
    for (size_t pTindex = 0; pTindex < pT.size(); pTindex++) {
        cout << pT[pTindex] << "\t";

        double total = 0;
        for (size_t hfgindex = 0; hfgindex < hfgroups.size(); hfgindex++) {
            rc->result(pTindex, hfgindex, &l_real, &l_imag, &l_error);
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
    }
    
    for (vector<HardFactorList*>::iterator it = hfgroups.begin(); it != hfgroups.end(); it++) {
        delete *it;
    }
    delete gctx.gdist;
    delete cpl;
    return 0;
}
