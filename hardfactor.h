#ifndef _HARD_FACTOR_H_
#define _HARD_FACTOR_H_

#include "integrationcontext.h"
#include "integrationtype.h"

class HardFactorTerm;

class HardFactor {
public:
    virtual const char* get_name() const = 0;
    virtual const size_t get_term_count() const = 0;
    virtual const HardFactorTerm* const* get_terms() const = 0;
};

class HardFactorTerm : public HardFactor {
public:
    HardFactorTerm() : p_this(this) {};
    virtual const IntegrationType* get_type() const = 0;
    virtual void Fs(const IntegrationContext* ictx, double* real, double* imag) const { *real = 0; *imag = 0; }
    virtual void Fn(const IntegrationContext* ictx, double* real, double* imag) const { *real = 0; *imag = 0; }
    virtual void Fd(const IntegrationContext* ictx, double* real, double* imag) const { *real = 0; *imag = 0; }
    const size_t get_term_count() const {
        return 1;
    }
    const HardFactorTerm* const* get_terms() const {
        return &p_this;
    }
private:
    const HardFactorTerm* p_this; // wtf why do I have to do this
};

#endif // _HARD_FACTOR_H_
