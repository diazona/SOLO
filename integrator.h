#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include <map>
#include <vector>
#include "context.h"
#include "integrationcontext.h"
#include "integrationtype.h"
#include "hardfactor.h"

typedef enum {MC_PLAIN, MC_MISER, MC_VEGAS} integration_strategy;

typedef std::vector<const HardFactor*> HardFactorList;
typedef std::vector<const HardFactorTerm*> HardFactorTermList;
typedef std::map<const IntegrationType*, HardFactorTermList> HardFactorTypeMap;

class Integrator {
private:
    IntegrationType const* current_type;
    IntegrationContext ictx;
    HardFactorTypeMap terms;
    integration_strategy strategy;
    void (*callback)(const IntegrationContext*, double, double);
public:
    Integrator(const Context* ctx, const ThreadLocalContext* tlctx, integration_strategy strategy, HardFactorList hflist);
    ~Integrator();
    void update1D(const double* values);
    void update2D(const double* values);
    void evaluate_1D_integrand(double* real, double* imag);
    void evaluate_2D_integrand(double* real, double* imag);
    void integrate(double* real, double* imag, double* error);
    void set_current_integration_type(IntegrationType const* new_type) {
        current_type = new_type;
    }
    void set_callback(void (*callback)(const IntegrationContext*, double, double)) {
        this->callback = callback;
    }
private:
    void integrate_impl(const size_t core_dimensions, double* result, double* error);
};

#endif // _INTEGRATOR_H_
