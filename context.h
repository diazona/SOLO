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

    Context(double x0, double A, double c, double lambda, double mu2, double Nc, double Nf, double CF, double TR, double Sperp, double pT2, double sqs, double Y, GluonDistribution* gdist, Coupling* cpl) :
     x0(x0), A(A), c(c), lambda(lambda), mu2(mu2), Nc(Nc), Nf(Nf), CF(CF), TR(TR), Sperp(Sperp), pT2(pT2), sqs(sqs), Y(Y), gdist(gdist), cpl(cpl) {
         recalculate();
    }

    void recalculate() {
         Q02x0lambda = c * pow(A, 1.0d/3.0d) * pow(x0, lambda);
         tau = sqrt(pT2)/sqs*exp(Y);
    }
};

class ThreadLocalContext {
    friend class IntegrationContext;
    c_mstwpdf* pdf_object;
    DSSpiNLO* ff_object;

public:
    ThreadLocalContext(const char* pdf_filename, const char* ff_filename) : pdf_object(new c_mstwpdf(pdf_filename)), ff_object(new DSSpiNLO(ff_filename)) {}
    
    ~ThreadLocalContext() {
        if (pdf_object) {
            delete pdf_object;
        }
        if (ff_object) {
            delete ff_object;
        }
    }
};

#endif // _CONTEXT_H_
