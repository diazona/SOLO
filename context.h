#ifndef _CONTEXT_H_
#define _CONTEXT_H_

#include <iterator>
#include <string>
#include <vector>
#include "mstwpdf.h"
#include "dss_pinlo.h"
#include "coupling.h"
#include "gluondist.h"

class Context {
public:
    const double x0;
    const double A;
    const double c;
    const double lambda;
    const double lambdaQCD;
    const double mu2;
    const double Nc;
    const double Nf;
    const double CF;
    const double TR;
    const double Sperp;
    const double pT2;
    const double sqs;
    const double Y;
    const std::string pdf_filename;
    const std::string ff_filename;

    GluonDistribution* gdist;
    Coupling* cpl;

    const double Q02x0lambda;
    const double tau;

    Context(
        const double x0,
        const double A,
        const double c,
        const double lambda,
        const double lambdaQCD,
        const double mu2,
        const double Nc,
        const double Nf,
        const double CF,
        const double TR,
        const double Sperp,
        const double pT2,
        const double sqs,
        const double Y,
        const std::string pdf_filename,
        const std::string ff_filename,
        GluonDistribution* gdist,
        Coupling* cpl) :
     x0(x0),
     A(A),
     c(c),
     lambda(lambda),
     lambdaQCD(lambdaQCD),
     mu2(mu2),
     Nc(Nc),
     Nf(Nf),
     CF(CF),
     TR(TR),
     Sperp(Sperp),
     pT2(pT2),
     sqs(sqs),
     Y(Y),
     pdf_filename(pdf_filename),
     ff_filename(ff_filename),
     gdist(gdist),
     cpl(cpl),
     Q02x0lambda(c * pow(A, 1.0d/3.0d) * pow(x0, lambda)),
     tau(sqrt(pT2)/sqs*exp(Y)) {}
};

std::ostream operator<<(std::ostream out, Context& ctx);

class ContextCollectionIterator;

class ContextCollection {
public:
    static const double unset = -1;
    typedef ContextCollectionIterator iterator;
    ContextCollection() :
      x0(0.000304),
      A(unset),
      c(unset),
      lambda(0.288),
      lambdaQCD(unset),
      mu2(unset),
      Nc(3),
      Nf(3),
      CF(1.5),
      TR(0.5),
      Sperp(unset),
      sqs(unset),
      pdf_filename("mstw2008nlo.00.dat"),
      ff_filename("PINLO.DAT"),
      gdist(NULL),
      cpl(NULL) {
    }
    ContextCollection(const std::string& filename) :
      x0(0.000304),
      A(unset),
      c(unset),
      lambda(0.288),
      lambdaQCD(unset),
      mu2(unset),
      Nc(3),
      Nf(3),
      CF(1.5),
      TR(0.5),
      Sperp(unset),
      sqs(unset),
      pdf_filename("mstw2008nlo.00.dat"),
      ff_filename("PINLO.DAT"),
      gdist(NULL),
      cpl(NULL) {
        ifstream in(filename.c_str());
        read_config(in);
        in.close();
    }
    ~ContextCollection() {
        delete gdist;
        gdist = NULL;
        delete cpl;
        cpl = NULL;
    }
    Context get_context(size_t n);
    Context operator[](size_t n) {
        return get_context(n);
    }
    size_t size() {
        return pT2.size() * Y.size();
    }
    iterator begin();
    iterator end();
    
    friend std::istream& operator>>(std::istream& in, ContextCollection& cc);
    friend std::ostream& operator<<(std::ostream& out, ContextCollection& cc);
    friend class ThreadLocalContext;
    
    double x0;
    double A;
    double c;
    double lambda;
    double lambdaQCD;
    double mu2;
    double Nc;
    double Nf;
    double CF;
    double TR;
    double Sperp;
    std::vector<double> pT2;
    double sqs;
    std::vector<double> Y;
    std::string pdf_filename;
    std::string ff_filename;
    GluonDistribution* gdist;
    Coupling* cpl;
    
    double Q02x0lambda() {
        return c * pow(A, 1.0d/3.0d) * pow(x0, lambda);
    }
private:
    void read_config(std::istream& in);
};

class ContextCollectionIterator : public iterator<random_access_iterator_tag, Context> {
    ContextCollection& cc;
    size_t n;
public:
    ContextCollectionIterator(ContextCollection& cc, size_t n) : cc(cc), n(n) {}
    ContextCollectionIterator(const ContextCollectionIterator& other) : cc(other.cc), n(other.n) {}
    ContextCollectionIterator& operator++() {
        ++n;
        return *this;
    }
    ContextCollectionIterator operator++(int) {
        ContextCollectionIterator tmp(*this);
        operator++();
        return tmp;
    }
    bool operator==(const ContextCollectionIterator& rhs) {
        return n == rhs.n && &cc == &(rhs.cc);
    }
    bool operator!=(const ContextCollectionIterator& rhs) {
        return n != rhs.n || &cc != &(rhs.cc);
    }
    Context operator*() {
        return cc.get_context(n);
    }
};

std::istream& operator>>(std::istream& in, ContextCollection& cc);
std::ostream& operator<<(std::ostream& out, ContextCollection& cc);

class ThreadLocalContext {
    friend class IntegrationContext;
    c_mstwpdf* pdf_object;
    DSSpiNLO* ff_object;

public:
    ThreadLocalContext(const char* pdf_filename, const char* ff_filename) :
      pdf_object(new c_mstwpdf(pdf_filename)),
      ff_object(new DSSpiNLO(ff_filename)) {};
    ThreadLocalContext(const Context& ctx) :
      pdf_object(new c_mstwpdf(ctx.pdf_filename.c_str())),
      ff_object(new DSSpiNLO(ctx.ff_filename.c_str())) {};
    ThreadLocalContext(const ContextCollection& ctx) :
      pdf_object(new c_mstwpdf(ctx.pdf_filename.c_str())),
      ff_object(new DSSpiNLO(ctx.ff_filename.c_str())) {};
    
    ~ThreadLocalContext() {
        delete pdf_object;
        delete ff_object;
    }
};

#endif // _CONTEXT_H_
