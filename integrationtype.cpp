#include "integrationtype.h"

static const double inf = 10;

void PlainIntegrationType::fill_min(IntegrationContext& ictx, const size_t core_dimensions, double* min) const {
    size_t i = 0;
    while (i < core_dimensions) { min[i++] = ictx.ctx->tau; }
    while (i < core_dimensions + extra_dimensions) { min[i++] = -inf; }
}
void PlainIntegrationType::fill_max(IntegrationContext& ictx, const size_t core_dimensions, double* max) const {
    size_t i = 0;
    while (i < core_dimensions) { max[i++] = 1; }
    while (i < core_dimensions + extra_dimensions) { max[i++] = inf; }
}

void DipoleIntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    if (core_dimensions == 1) {
        ictx.update_parton_factors(values[0], 1);
        ictx.update_positions(values[1], values[2], 0, 0, 0, 0);
    }
    else {
        ictx.update_parton_factors(values[0], values[1]);
        ictx.update_positions(values[2], values[3], 0, 0, 0, 0);
    }
}

void QuadrupoleIntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    if (core_dimensions == 1) {
        ictx.update_parton_factors(values[0], 1);
        ictx.update_positions(values[1], values[2], values[3], values[4], 0, 0);
    }
    else {
        ictx.update_parton_factors(values[0], values[1]);
        ictx.update_positions(values[2], values[3], values[4], values[5], 0, 0);
    }
}

void Momentum1IntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    if (core_dimensions == 1) {
        ictx.update_parton_factors(values[0], 1);
        ictx.update_momenta(values[1], values[2], 0, 0, 0, 0);
    }
    else {
        ictx.update_parton_factors(values[0], values[1]);
        ictx.update_momenta(values[2], values[3], 0, 0, 0, 0);
    }
}

void Momentum2IntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    if (core_dimensions == 1) {
        ictx.update_parton_factors(values[0], 1);
        ictx.update_momenta(values[1], values[2], values[3], values[4], 0, 0);
    }
    else {
        ictx.update_parton_factors(values[0], values[1]);
        ictx.update_momenta(values[2], values[3], values[4], values[5], 0, 0);
    }
}

void Momentum3IntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    if (core_dimensions == 1) {
        ictx.update_parton_factors(values[0], 1);
        ictx.update_momenta(values[1], values[2], values[3], values[4], values[5], values[6]);
    }
    else {
        ictx.update_parton_factors(values[0], values[1]);
        ictx.update_momenta(values[2], values[3], values[4], values[5], values[6], values[7]);
    }
}

void XiPIntegrationType::fill_min(IntegrationContext& ictx, const size_t core_dimensions, double* min) const {
    size_t i = 0;
    while (i < core_dimensions) { min[i++] = ictx.ctx->tau; }
    min[i++] = 0;
    while (i < core_dimensions + extra_dimensions) { min[i++] = -inf; }
}
void XiPIntegrationType::fill_max(IntegrationContext& ictx, const size_t core_dimensions, double* max) const {
    size_t i = 0;
    while (i < core_dimensions) { max[i++] = 1; }
    max[i++] = 1;
    while (i < core_dimensions + extra_dimensions) { max[i++] = inf; }
}

void Momentum1XiPIntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    if (core_dimensions == 1) {
        ictx.update_parton_factors(values[0], 1);
        ictx.update_auxiliary(values[1]);
        ictx.update_momenta(values[2], values[3], 0, 0, 0, 0);
    }
    else {
        ictx.update_parton_factors(values[0], values[1]);
        ictx.update_auxiliary(values[2]);
        ictx.update_momenta(values[3], values[4], 0, 0, 0, 0);
    }
}

void Momentum2XiPIntegrationType::update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const {
    if (core_dimensions == 1) {
        ictx.update_parton_factors(values[0], 1);
        ictx.update_auxiliary(values[1]);
        ictx.update_momenta(values[2], values[3], values[4], values[5], 0, 0);
    }
    else {
        ictx.update_parton_factors(values[0], values[1]);
        ictx.update_auxiliary(values[2]);
        ictx.update_momenta(values[3], values[4], values[5], values[6], 0, 0);
    }
}
