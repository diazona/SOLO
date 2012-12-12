#ifndef _INTEGRATOR_H_
#define _INTEGRATOR_H_

#include <map>
#include <vector>
#include "context.h"
#include "integrationcontext.h"
#include "integrationtype.h"
#include "hardfactor.h"

typedef enum {MC_PLAIN, MC_MISER, MC_VEGAS} integration_strategy;

typedef std::vector<HardFactor*> HardFactorList;
typedef std::map<IntegrationType*, HardFactorList> HardFactorTypeMap;

class Integrator {
private:
    IntegrationType* current_type;
    IntegrationContext ictx;
    HardFactorTypeMap terms;
    integration_strategy strategy;
    void (*callback)(IntegrationContext*, double, double);
public:
    Integrator(Context* ctx, integration_strategy strategy, HardFactorList hflist);
    ~Integrator();
    void update1D(double* values);
    void update2D(double* values);
    void evaluate_1D_integrand(double* real, double* imag);
    void evaluate_2D_integrand(double* real, double* imag);
    void integrate(double* real, double* imag, double* error);
    void set_current_integration_type(IntegrationType* new_type) {
        current_type = new_type;
    }
    void set_callback(void (*callback)(IntegrationContext*, double, double)) {
        this->callback = callback;
    }
private:
    void integrate_impl(size_t core_dimensions, double* result, double* error);
};

#endif // _INTEGRATOR_H_
