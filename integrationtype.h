#ifndef _INTEGRATIONTYPE_H_
#define _INTEGRATIONTYPE_H_

#include "integrationcontext.h"

class IntegrationType {
public:
    // extra_dimensions is the number of integrals NOT including z and y
    const size_t extra_dimensions;
    IntegrationType(const size_t extra_dimensions) : extra_dimensions(extra_dimensions) {}
    virtual void fill_min(IntegrationContext& ictx, const size_t core_dimensions, double* min) const = 0;
    virtual void fill_max(IntegrationContext& ictx, const size_t core_dimensions, double* max) const = 0;
    virtual void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const = 0;
};

class PlainIntegrationType : public IntegrationType {
protected:
    PlainIntegrationType(const size_t extra_dimensions) : IntegrationType(extra_dimensions) {}
    void fill_min(IntegrationContext& ictx, const size_t core_dimensions, double* min) const;
    void fill_max(IntegrationContext& ictx, const size_t core_dimensions, double* max) const;
};

class DipoleIntegrationType : public PlainIntegrationType {
public:
    static const DipoleIntegrationType* get_instance() {
        static DipoleIntegrationType instance;
        return &instance;
    }
private:
    DipoleIntegrationType() : PlainIntegrationType(2) {}
    DipoleIntegrationType(DipoleIntegrationType const&);
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

class QuadrupoleIntegrationType : public PlainIntegrationType {
public:
    static const QuadrupoleIntegrationType* get_instance() {
        static QuadrupoleIntegrationType instance;
        return &instance;
    }
private:
    QuadrupoleIntegrationType() : PlainIntegrationType(4) {}
    QuadrupoleIntegrationType(QuadrupoleIntegrationType const&);
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

class Momentum1IntegrationType : public PlainIntegrationType {
public:
    static const Momentum1IntegrationType* get_instance() {
        static Momentum1IntegrationType instance;
        return &instance;
    }
private:
    Momentum1IntegrationType() : PlainIntegrationType(2) {}
    Momentum1IntegrationType(Momentum1IntegrationType const&);
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

class Momentum2IntegrationType : public PlainIntegrationType {
public:
    static const Momentum2IntegrationType* get_instance() {
        static Momentum2IntegrationType instance;
        return &instance;
    }
private:
    Momentum2IntegrationType() : PlainIntegrationType(4) {}
    Momentum2IntegrationType(Momentum2IntegrationType const&);
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

class Momentum3IntegrationType : public PlainIntegrationType {
public:
    static const Momentum3IntegrationType* get_instance() {
        static Momentum3IntegrationType instance;
        return &instance;
    }
private:
    Momentum3IntegrationType() : PlainIntegrationType(6) {}
    Momentum3IntegrationType(Momentum3IntegrationType const&);
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

class XiPIntegrationType : public IntegrationType {
protected:
    XiPIntegrationType(size_t extra_dimensions) : IntegrationType(extra_dimensions) {}
    void fill_min(IntegrationContext& ictx, const size_t core_dimensions, double* min) const;
    void fill_max(IntegrationContext& ictx, const size_t core_dimensions, double* max) const;
};

class Momentum1XiPIntegrationType : public XiPIntegrationType {
public:
    static const Momentum1XiPIntegrationType* get_instance() {
        static Momentum1XiPIntegrationType instance;
        return &instance;
    }
private:
    Momentum1XiPIntegrationType() : XiPIntegrationType(3) {}
    Momentum1XiPIntegrationType(Momentum1XiPIntegrationType const&);
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

class Momentum2XiPIntegrationType : public XiPIntegrationType {
public:
    static const Momentum2XiPIntegrationType* get_instance() {
        static Momentum2XiPIntegrationType instance;
        return &instance;
    }
private:
    Momentum2XiPIntegrationType() : XiPIntegrationType(5) {}
    Momentum2XiPIntegrationType(Momentum2XiPIntegrationType const&);
    void update(IntegrationContext& ictx, const size_t core_dimensions, const double* values) const;
};

extern const double inf;

#endif
