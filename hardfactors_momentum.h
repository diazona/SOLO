#ifndef _HARDFACTORS_MOMENTUM_H_
#define _HARDFACTORS_MOMENTUM_H_

// Hard factors in momentum space

#include "hardfactor.h"

namespace momentum {

// class H02qq : public HardFactor {
// public:
//     const char* get_name() {
//         return "H02qq";
//     }
//     IntegrationType* get_type() {
//         return dipole;
//     }
//     void Fd(IntegrationContext* ictx, double* real, double* imag);
// };
// 
// class H12qq : public HardFactor {
// public:
//     const char* get_name() {
//         return "H12qq";
//     }
//     IntegrationType* get_type() {
//         return dipole;
//     }
//     void Fs(IntegrationContext* ictx, double* real, double* imag);
//     void Fd(IntegrationContext* ictx, double* real, double* imag);
// };
// 
// class H14qq : public HardFactor {
// public:
//     const char* get_name() {
//         return "H14qq";
//     }
//     IntegrationType* get_type() {
//         return quadrupole;
//     }
//     void Fs(IntegrationContext* ictx, double* real, double* imag);
//     void Fd(IntegrationContext* ictx, double* real, double* imag);
// };
// 
// class H14qqResidual : public HardFactor {
// public:
//     const char* get_name() {
//         return "H14qqResidual";
//     }
//     IntegrationType* get_type() {
//         return dipole;
//     }
//     void Fd(IntegrationContext* ictx, double* real, double* imag);
// };
// 
// class H02gg : public HardFactor {
// public:
//     const char* get_name() {
//         return "H02gg";
//     }
//     IntegrationType* get_type() {
//         return dipole;
//     }
//     void Fd(IntegrationContext* ictx, double* real, double* imag);
// };
// 
// class H12gg : public HardFactor {
// public:
//     const char* get_name() {
//         return "H12gg";
//     }
//     IntegrationType* get_type() {
//         return dipole;
//     }
//     void Fs(IntegrationContext* ictx, double* real, double* imag);
//     void Fn(IntegrationContext* ictx, double* real, double* imag);
//     void Fd(IntegrationContext* ictx, double* real, double* imag);
// };

class H12qqbar : public HardFactor {
public:
    const char* get_name() {
        return "H12qqbar";
    }
    IntegrationType* get_type() {
        return Momentum1XiPIntegrationType::get_instance();
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag);
};

class H16ggSingular : public HardFactor {
public:
    const char* get_name() {
        return "H16gg";
    }
    IntegrationType* get_type() {
        return Momentum3IntegrationType::get_instance();
    }
    void Fs(IntegrationContext* ictx, double* real, double* imag);
};

class H16ggDelta : public HardFactor {
public:
    const char* get_name() {
        return "H16gg";
    }
    IntegrationType* get_type() {
        return Momentum2XiPIntegrationType::get_instance();
    }
    void Fd(IntegrationContext* ictx, double* real, double* imag);
};

// class H112gq : public HardFactor {
// public:
//     const char* get_name() {
//         return "H112gq";
//     }
//     IntegrationType* get_type() {
//         return dipole;
//     }
//     void Fn(IntegrationContext* ictx, double* real, double* imag);
// };
// 
// class H122gq : public HardFactor {
// public:
//     const char* get_name() {
//         return "H122gq";
//     }
//     IntegrationType* get_type() {
//         return dipole;
//     }
//     void Fn(IntegrationContext* ictx, double* real, double* imag);
// };

class H14gq : public HardFactor {
public:
    const char* get_name() {
        return "H14gq";
    }
    IntegrationType* get_type() {
        return Momentum2IntegrationType::get_instance();
    }
    void Fn(IntegrationContext* ictx, double* real, double* imag);
};

// class H112qg : public HardFactor {
// public:
//     const char* get_name() {
//         return "H112qg";
//     }
//     IntegrationType* get_type() {
//         return dipole;
//     }
//     void Fn(IntegrationContext* ictx, double* real, double* imag);
// };
// 
// class H122qg : public HardFactor {
// public:
//     const char* get_name() {
//         return "H122qg";
//     }
//     IntegrationType* get_type() {
//         return dipole;
//     }
//     void Fn(IntegrationContext* ictx, double* real, double* imag);
// };

class H14qg : public HardFactor {
public:
    const char* get_name() {
        return "H14qg";
    }
    IntegrationType* get_type() {
        return Momentum2IntegrationType::get_instance();
    }
    void Fn(IntegrationContext* ictx, double* real, double* imag);
};

}

#endif // _HARDFACTORS_MOMENTUM_H_