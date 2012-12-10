// Hard factors in position space

#include "hardfactor.h"

class H02qq : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag);
};

class H12qq : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fs(IntegrationContext* ictx, double* real, double* imag);
    void Fd(IntegrationContext* ictx, double* real, double* imag);
};

class H14qq : public HardFactor {
public:
    term_type get_type() {
        return quadrupole;
    }
    void Fs(IntegrationContext* ictx, double* real, double* imag);
    void Fd(IntegrationContext* ictx, double* real, double* imag);
};

class H14qqResidual : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag);
};

class H02gg : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag);
};

class H12gg : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fs(IntegrationContext* ictx, double* real, double* imag);
    void Fn(IntegrationContext* ictx, double* real, double* imag);
    void Fd(IntegrationContext* ictx, double* real, double* imag);
};

class H12qqbar : public HardFactor {
public:
    term_type get_type() {
        return quadrupole;
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag);
};

class H12qqbarResidual : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag);
};

class H16gg : public HardFactor {
public:
    term_type get_type() {
        return quadrupole;
    }
    void Fs(IntegrationContext* ictx, double* real, double* imag);
    void Fd(IntegrationContext* ictx, double* real, double* imag);
};

class H16ggResidual : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag);
};

class H112gq : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fn(IntegrationContext* ictx, double* real, double* imag);
};

class H122gq : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fn(IntegrationContext* ictx, double* real, double* imag);
};

class H14gq : public HardFactor {
public:
    term_type get_type() {
        return quadrupole;
    }
    void Fn(IntegrationContext* ictx, double* real, double* imag);
};

class H112qg : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fn(IntegrationContext* ictx, double* real, double* imag);
};

class H122qg : public HardFactor {
public:
    term_type get_type() {
        return dipole;
    }
    void Fn(IntegrationContext* ictx, double* real, double* imag);
};

class H14qg : public HardFactor {
public:
    term_type get_type() {
        return quadrupole;
    }
    void Fn(IntegrationContext* ictx, double* real, double* imag);
};

