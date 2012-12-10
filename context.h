#ifndef _CONTEXT_H_
#define _CONTEXT_H_

#include "mstwpdf.h"
#include "dss_pinlo.h"
#include "coupling.h"
#include "gluondist.h"

class Context {
public:
    double x0;
    double A;
    double c;
    double lambda;
    double mu2;
    double Nc;
    double Nf;
    double CF;
    double TR;
    double Sperp;
    double pT2;
    double sqs;
    double Y;
    double Q02x0lambda;
    double tau;
    
    GluonDistribution* gdist;
    Coupling* cpl;

    // NOTE: these contain state that should be associated with IntegrationContext.
    // So don't use one Context with more than one IntegrationContext at once.
    c_mstwpdf* pdf_object;
    DSSpiNLO* ff_object;

    Context(double x0, double A, double c, double lambda, double mu2, double Nc, double Nf, double CF, double TR, double Sperp, double pT2, double sqs, double Y, GluonDistribution* gdist, Coupling* cpl, const char* pdf_filename, const char* ff_filename) :
     x0(x0), A(A), c(c), lambda(lambda), mu2(mu2), Nc(Nc), Nf(Nf), CF(CF), TR(TR), Sperp(Sperp), pT2(pT2), sqs(sqs), Y(Y), gdist(gdist), cpl(cpl) {
         recalculate();
         pdf_object = new c_mstwpdf(pdf_filename);
         ff_object = new DSSpiNLO(ff_filename);
    }

    ~Context() {
        if (pdf_object) {
            delete pdf_object;
        }
        if (ff_object) {
            delete ff_object;
        }
    }

    void recalculate() {
         Q02x0lambda = c * pow(A, 1.0d/3.0d) * pow(x0, lambda);
         tau = sqrt(pT2)/sqs*exp(Y);
    }
};

#endif // _CONTEXT_H_
