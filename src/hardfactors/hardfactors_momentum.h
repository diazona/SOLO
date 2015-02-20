/*
 * Part of oneloopcalc
 * 
 * Copyright 2012 David Zaslavsky
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _HARDFACTORS_MOMENTUM_H_
#define _HARDFACTORS_MOMENTUM_H_

// Hard factors in momentum space

#include "hardfactor.h"

namespace momentum {

extern const NoIntegrationType none;
extern const MomentumIntegrationType momentum1;
extern const MomentumIntegrationType momentum2;
extern const MomentumIntegrationType momentum3;
extern const RadialMomentumIntegrationType radialmomentum1;
extern const RadialMomentumIntegrationType radialmomentum2;
extern const RadialMomentumIntegrationType radialmomentum3;
extern const XiPIntegrationType momentumxip1;
extern const XiPIntegrationType momentumxip2;
extern const QLimitedMomentumIntegrationType qlim;

class H02qq : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H02qq";
    }
    const IntegrationType* get_type() const {
        return &none;
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

// class H12qq : public HardFactorTerm {
// public:
//     const char* get_name() {
//         return "H12qq";
//     }
//     IntegrationType* get_type() {
//         return dipole;
//     }
//     void Fs(const IntegrationContext* ictx, double* real, double* imag);
//     void Fd(const IntegrationContext* ictx, double* real, double* imag);
// };

class H14qqSingular : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H14qqSingular";
    }
    const IntegrationType* get_type() const {
        return &momentum2;
    }
    void Fs(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H14qqDelta : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H14qqDelta";
    }
    const IntegrationType* get_type() const {
        return &momentumxip1;
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H14qq : public HardFactor {
public:
    H14qq() {
        terms[0] = new H14qqSingular();
        terms[1] = new H14qqDelta();
    }
    ~H14qq() {
        delete terms[0];
        terms[0] = NULL;
        delete terms[1];
        terms[1] = NULL;
    }
    const char* get_name() const {
        return "H14qq";
    }
    const size_t get_term_count() const {
        return 2;
    }
    const HardFactorTerm* const* get_terms() const {
        return terms;
    }
private:
    const HardFactorTerm* terms[2];
};

class H1qqCorrectionA : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H1qqCoA";
    }
    const IntegrationType* get_type() const {
        return &momentum1;
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H1qqCorrectionB : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H1qqCoB";
    }
    const IntegrationType* get_type() const {
        return &momentum1;
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H1qqCorrection : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H1qqCo";
    }
    const IntegrationType* get_type() const {
        return &momentum1;
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H1qqExact : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H1qqExact";
    }
    const IntegrationType* get_type() const {
        return &qlim;
    }
    void Fn(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H02gg : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H02gg";
    }
    const IntegrationType* get_type() const {
        return &momentum1;
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

// class H12gg : public HardFactorTerm {
// public:
//     const char* get_name() {
//         return "H12gg";
//     }
//     IntegrationType* get_type() {
//         return dipole;
//     }
//     void Fs(const IntegrationContext* ictx, double* real, double* imag);
//     void Fn(const IntegrationContext* ictx, double* real, double* imag);
//     void Fd(const IntegrationContext* ictx, double* real, double* imag);
// };

class H12qqbar : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H12qqbar";
    }
    const IntegrationType* get_type() const {
        return &momentumxip1;
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H16ggSingular : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H16ggSingular";
    }
    const IntegrationType* get_type() const {
        return &momentum3;
    }
    void Fs(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H16ggDelta : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H16ggDelta";
    }
    const IntegrationType* get_type() const {
        return &momentumxip2;
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H16gg : public HardFactor {
public:
    H16gg() {
        terms[0] = new H16ggSingular();
        terms[1] = new H16ggDelta();
    }
    ~H16gg() {
        delete terms[0];
        terms[0] = NULL;
        delete terms[1];
        terms[1] = NULL;
    }
    const char* get_name() const {
        return "H16gg";
    }
    const size_t get_term_count() const {
        return 2;
    }
    const HardFactorTerm* const* get_terms() const {
        return terms;
    }
private:
    const HardFactorTerm* terms[2];
};

class H1ggCorrection : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H1ggCo";
    }
    const IntegrationType* get_type() const {
        return &momentum1;
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H1ggExact : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H1ggExact";
    }
    const IntegrationType* get_type() const {
        return &qlim;
    }
    void Fn(const IntegrationContext* ictx, double* real, double* imag) const;
};


// class H112gq : public HardFactorTerm {
// public:
//     const char* get_name() {
//         return "H112gq";
//     }
//     IntegrationType* get_type() {
//         return dipole;
//     }
//     void Fn(const IntegrationContext* ictx, double* real, double* imag);
// };
// 
// class H122gq : public HardFactorTerm {
// public:
//     const char* get_name() {
//         return "H122gq";
//     }
//     IntegrationType* get_type() {
//         return dipole;
//     }
//     void Fn(const IntegrationContext* ictx, double* real, double* imag);
// };

class H14gq : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H14gq";
    }
    const IntegrationType* get_type() const {
        return &momentum2;
    }
    void Fn(const IntegrationContext* ictx, double* real, double* imag) const;
};

// class H112qg : public HardFactorTerm {
// public:
//     const char* get_name() {
//         return "H112qg";
//     }
//     IntegrationType* get_type() {
//         return dipole;
//     }
//     void Fn(const IntegrationContext* ictx, double* real, double* imag);
// };
// 
// class H122qg : public HardFactorTerm {
// public:
//     const char* get_name() {
//         return "H122qg";
//     }
//     IntegrationType* get_type() {
//         return dipole;
//     }
//     void Fn(const IntegrationContext* ictx, double* real, double* imag);
// };

class H14qg : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H14qg";
    }
    const IntegrationType* get_type() const {
        return &momentum2;
    }
    void Fn(const IntegrationContext* ictx, double* real, double* imag) const;
};

class registry : public HardFactorRegistry {
public:
    static registry* get_instance() {
        static registry instance;
        return &instance;
    }
private:
    registry() {}
};

}

#endif // _HARDFACTORS_MOMENTUM_H_
