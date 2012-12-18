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
#include <iostream>
#include <cstdlib>
#include <cstring>
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

int main(int argc, char** argv) {
    bool separate = false;
    GluonDistribution* gdist = new GBWGluonDistribution();
    Coupling* cpl = new FixedCoupling(0.2 / (2*M_PI));
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
    }
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
      gdist,
      cpl,
      "mstw2008nlo.00.dat", "PINLO.DAT");

    double pT[] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4};
    vector<double> pTlist(pT, pT + sizeof(pT)/sizeof(double));
    size_t pTlen = sizeof(pT)/sizeof(double);
    position::H02qq         _h02qq;
    position::H02gg         _h02gg;
    position::H12qq         _h12qq;
    position::H14qq         _h14qq;
    position::H12gg         _h12gg;
    momentum::H12qqbar      _h12qqbar;
//     position::H12qqbar      _h12qqbar;
    momentum::H16gg         _h16gg;
//     position::H16gg         _h16gg;
    position::H112gq        _h112gq;
    position::H122gq        _h122gq;
    momentum::H14gq         _h14gq;
//     position::H14gq         _h14gq;
    position::H112qg        _h112qg;
    position::H122qg        _h122qg;
    momentum::H14qg         _h14qg;
//     position::H14qg         _h14qg;
    
    vector<HardFactorList*> hfgroups;
    if (separate) {
        hfgroups.push_back(new HardFactorList(1, &_h02qq));
        hfgroups.push_back(new HardFactorList(1, &_h02gg));
        hfgroups.push_back(new HardFactorList(1, &_h12qq));
        hfgroups.push_back(new HardFactorList(1, &_h14qq));
        hfgroups.push_back(new HardFactorList(1, &_h12gg));
        hfgroups.push_back(new HardFactorList(1, &_h12qqbar));
        hfgroups.push_back(new HardFactorList(1, &_h16gg));
        hfgroups.push_back(new HardFactorList(1, &_h112gq));
        hfgroups.push_back(new HardFactorList(1, &_h122gq));
        hfgroups.push_back(new HardFactorList(1, &_h14gq));
        hfgroups.push_back(new HardFactorList(1, &_h112qg));
        hfgroups.push_back(new HardFactorList(1, &_h122qg));
        hfgroups.push_back(new HardFactorList(1, &_h14qg));
    }
    else {
        HardFactorList* lo = new HardFactorList();
            lo->push_back(&_h02qq);
            lo->push_back(&_h02gg);
        hfgroups.push_back(lo);
        HardFactorList* nlo = new HardFactorList();
            nlo->push_back(&_h12qq);
            nlo->push_back(&_h14qq);
            nlo->push_back(&_h12gg);
            nlo->push_back(&_h12qqbar);
            nlo->push_back(&_h16gg);
            nlo->push_back(&_h112gq);
            nlo->push_back(&_h122gq);
            nlo->push_back(&_h14gq);
            nlo->push_back(&_h112qg);
            nlo->push_back(&_h122qg);
            nlo->push_back(&_h14qg);
        hfgroups.push_back(nlo);
    }
    
    ResultsCalculator* rc = new ResultsCalculator(&gctx, pTlist, hfgroups);
    rc->calculate();
    if (separate) {
        cout << "pT\th02qq\th02gg\th12qq\th14qq\th12gg\th12qqbar\th16gg\th112gq\th122gq\th14gq\th112qg\th122qg\th14qg\tlo\tnlo\tlo+nlo" << endl;
    }
    else {
        cout << "pT\tlo\tnlo\tlo+nlo" << endl;
    }
    double l_real, l_imag, l_error;
    for (size_t pTindex = 0; pTindex < pTlist.size(); pTindex++) {
        if (separate) {
            cout << pT[pTindex] << "\t";

            double lo = 0, nlo = 0;
            size_t hfgindex;
            for (hfgindex = 0; hfgindex < 2; hfgindex++) {
                rc->result(pTindex, hfgindex, &l_real, &l_imag, &l_error);
                cout << l_real << "±" << l_error << "\t";
                lo += l_real;
            }
            for (; hfgindex < hfgroups.size(); hfgindex++) {
                rc->result(pTindex, hfgindex, &l_real, &l_imag, &l_error);
                cout << l_real << "±" << l_error << "\t";
                nlo += l_real;
            }
            cout << lo << "\t" << nlo << "\t" << lo+nlo << endl;
        }
        else {
            double total = 0;
            cout << pT[pTindex] << "\t";
            rc->result(pTindex, 0, &l_real, &l_imag, &l_error);
            cout << l_real << "±" << l_error << "\t";
            total += l_real;
            rc->result(pTindex, 1, &l_real, &l_imag, &l_error);
            cout << l_real << "±" << l_error << "\t";
            total += l_real;
            cout << total << endl;
        }
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
    delete gdist;
    delete cpl;
    return 0;
}
