#ifndef _HARD_FACTOR_H_
#define _HARD_FACTOR_H_

#include "integrationcontext.h"

typedef enum {NONE=0, dipole=2, quadrupole=4} term_type;

class HardFactor {
public:
    virtual term_type get_type() = 0;
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