#ifndef _HARD_FACTOR_H_
#define _HARD_FACTOR_H_

#include "integrationcontext.h"
#include "integrationtype.h"

class HardFactor {
public:
    virtual const char* get_name() = 0;
    virtual IntegrationType* get_type() = 0;
    virtual void Fs(IntegrationContext* ictx, double* real, double* imag) { *real = 0; *imag = 0; }
    virtual void Fn(IntegrationContext* ictx, double* real, double* imag) { *real = 0; *imag = 0; }
    virtual void Fd(IntegrationContext* ictx, double* real, double* imag) { *real = 0; *imag = 0; }
};

class HardFactorGroup {
public:
    size_t hflen;
    HardFactor** hflist;
    HardFactorGroup(size_t hflen, HardFactor** hflist) : hflen(hflen), hflist(hflist) {};
};

#endif // _HARD_FACTOR_H_