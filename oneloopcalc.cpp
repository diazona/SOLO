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
#include "integrator.h"

using namespace std;

const int SUCCESS = 0;

static bool trace = false;

#ifdef DATAGRID
void write_hf_values(Context* ctx) {
    double real, imag;
    GBWGluonDistribution* gdist = new GBWGluonDistribution();
    IntegrationContext* ictx = new IntegrationContext(ctx);
    H02qq* hf02qq = new H02qq();
    H12qq* hf12qq = new H12qq();
    for (double z = ctx->tau; z <= 1; z += (1 - ctx->tau) / 20) {
        for (double rx = -5; rx <= 5; rx += 20./270.) {
            for (double ry = -5; ry <= 5; ry += 20./120.) {
                ictx->update(z, 1, rx, ry, 0, 0, 0, 0);
                cout << z << "\t" << rx << "\t" << ry << "\t"
                     << hf02qq->real_delta_contribution(ictx) << "\t" << hf02qq->imag_delta_contribution(ictx) << "\t"
                     << hf12qq->real_delta_contribution(ictx) << "\t" << hf12qq->imag_delta_contribution(ictx) << endl;
            }
        }
    }
}

int main(int argc, char** argv) {
    Context gctx(
      0.000304, // x0
      197,      // A
      0.56,     // c
      0.288,    // lambda
      10,       // mu2
      3,        // Nc
      3,        // Nf
      1.5,      // CF
      1.0,      // Sperp
      0.4,      // pT2
      200,      // sqs
      3.2,      // Y
      0.2 / (2*M_PI), // alphasbar
      "mstw2008nlo.00.dat", "PINLO.DAT");
    write_hf_values(&gctx);
    return 0;
}
#else
void write_data_point(IntegrationContext* ictx, double real, double imag) {
    if (ictx) {
//         if ((++count) % 500 == 0) {
            cout
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
        cout << endl;
    }
}

class ResultsCalculator {
private:
    Context* ctx;
    size_t pTlen;
    double* pTlist;
    size_t hfglen;
    HardFactorGroup** hfgroups;
    double* real;
    double* imag;
    double* error;
public:
    ResultsCalculator(Context* ctx, size_t pTlen, double* pTlist, size_t hfglen, HardFactorGroup** hfgroups) : ctx(ctx), pTlen(pTlen), pTlist(pTlist), hfglen(hfglen), hfgroups(hfgroups) {
        real = new double[hfglen * pTlen];
        imag = new double[hfglen * pTlen];
        error = new double[hfglen * pTlen];
    };
    ~ResultsCalculator() {
        delete[] real;
        delete[] imag;
        delete[] error;
    };
    size_t index_from(size_t pTindex, size_t hfgindex) {
        return pTindex * hfglen + hfgindex;
    };
    void result(size_t pTindex, size_t hfgindex, double* real, double* imag, double* error);
    void calculate();
};

void ResultsCalculator::result(size_t pTindex, size_t hfgindex, double* real, double* imag, double* error) {
    size_t index = index_from(pTindex, hfgindex);
    (*real) = this->real[index];
    (*imag) = this->imag[index];
    (*error) = this->error[index];
}

void ResultsCalculator::calculate() {
    for (size_t pTindex = 0; pTindex < pTlen; pTindex++) {
        ctx->pT2 = pTlist[pTindex] * pTlist[pTindex];
        ctx->recalculate();
        cerr << "Beginning calculation at pT = " << pTlist[pTindex] << endl;
        for (size_t hfgindex = 0; hfgindex < hfglen; hfgindex++) {
            size_t index = index_from(pTindex, hfgindex);
            Integrator* integrator = new Integrator(ctx, hfgroups[hfgindex]->hflen, hfgroups[hfgindex]->hflist);
            if (trace) {
                integrator->set_callback(write_data_point);
            }
            integrator->integrate(real + index, imag + index, error + index);
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
            integration_strategy = MC_MISER;
        }
        else if (strcmp(argv[i], "--vegas")==0) {
            integration_strategy = MC_VEGAS;
        }
        else if (strcmp(argv[i], "--trace")==0) {
            trace = true;
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
    size_t pTlen = sizeof(pT)/sizeof(double);
    HardFactor* hflist[] = {
        new H02qq(), // 0
        new H02gg(), // 1
        new H12qq(), // 2
        new H14qq(), new H14qqResidual(), // 3,4
        new H12gg(), // 5
        new H12qqbar(), new H12qqbarResidual(), // 6,7
        new H16gg(), new H16ggResidual(), // 8,9
        new H112gq(), // 10
        new H122gq(), // 11
        new H14gq(), // 12
        new H112qg(), // 13
        new H122qg(), // 14
        new H14qg()}; // 15
    HardFactorGroup* hfglist_separate[] = {
        new HardFactorGroup(1, hflist + 0),
        new HardFactorGroup(1, hflist + 1),
        new HardFactorGroup(1, hflist + 2),
        new HardFactorGroup(2, hflist + 3),
        new HardFactorGroup(1, hflist + 5),
        new HardFactorGroup(2, hflist + 6),
        new HardFactorGroup(2, hflist + 8),
        new HardFactorGroup(1, hflist + 10),
        new HardFactorGroup(1, hflist + 11),
        new HardFactorGroup(1, hflist + 12),
        new HardFactorGroup(1, hflist + 13),
        new HardFactorGroup(1, hflist + 14),
        new HardFactorGroup(1, hflist + 15)
    };
    HardFactorGroup* hfglist_total[] = {
        new HardFactorGroup(2, hflist + 0), // LO
        new HardFactorGroup(14, hflist + 2) // NLO
    };
    HardFactorGroup** hfglist = NULL;
    size_t hfglen;
    if (separate) {
        hfglist = hfglist_separate;
        hfglen = 13;
    }
    else {
        hfglist = hfglist_total;
        hfglen = 2;
    }
    
    ResultsCalculator* rc = new ResultsCalculator(&gctx, pTlen, pT, hfglen, hfglist);
    rc->calculate();
    if (separate) {
        cout << "pT\th02qq\th02gg\th12qq\th14qq\th12gg\th12qqbar\th16gg\th112gq\th122gq\th14gq\th112qg\th122qg\th14qg\tlo\tnlo\tlo+nlo" << endl;
    }
    else {
        cout << "pT\tlo\tnlo\tlo+nlo" << endl;
    }
    double l_real, l_imag, l_error;
    for (size_t pTindex = 0; pTindex < pTlen; pTindex++) {
        if (separate) {
            cout << pT[pTindex] << "\t";

            double lo = 0, nlo = 0;
            size_t hfgindex;
            for (hfgindex = 0; hfgindex < 2; hfgindex++) {
                rc->result(pTindex, hfgindex, &l_real, &l_imag, &l_error);
                cout << l_real << "±" << l_error << "\t";
                lo += l_real;
            }
            for (; hfgindex < hfglen; hfgindex++) {
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
    delete gdist;
    delete cpl;
    return 0;
}
#endif
