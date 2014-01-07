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

class H02qq : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H02qq";
    }
    const IntegrationType* get_type() const {
        return NoIntegrationType::get_instance();
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
        return Momentum2IntegrationType::get_instance();
    }
    void Fs(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H14qqDelta : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H14qqDelta";
    }
    const IntegrationType* get_type() const {
        return Momentum1XiPIntegrationType::get_instance();
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
        return Momentum1IntegrationType::get_instance();
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H1qqCorrectionB : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H1qqCoB";
    }
    const IntegrationType* get_type() const {
        return Momentum1IntegrationType::get_instance();
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H1qqCorrection : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H1qqCo";
    }
    const IntegrationType* get_type() const {
        return Momentum1IntegrationType::get_instance();
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H02gg : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H02gg";
    }
    const IntegrationType* get_type() const {
        return Momentum1IntegrationType::get_instance();
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
        return Momentum1XiPIntegrationType::get_instance();
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H16ggSingular : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H16ggSingular";
    }
    const IntegrationType* get_type() const {
        return Momentum3IntegrationType::get_instance();
    }
    void Fs(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H16ggDelta : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H16ggDelta";
    }
    const IntegrationType* get_type() const {
        return Momentum2XiPIntegrationType::get_instance();
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
        return Momentum1IntegrationType::get_instance();
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
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
        return Momentum2IntegrationType::get_instance();
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
        return Momentum2IntegrationType::get_instance();
    }
    void Fn(const IntegrationContext* ictx, double* real, double* imag) const;
};

class registry : public HardFactorRegistry {
public:
    static const registry* get_instance() {
        static registry instance;
        return &instance;
    }
private:
    registry() {
        add_hard_factor(new H02qq(), true);
//         add_hard_factor(new H12qq(), true);
        add_hard_factor(new H14qq(), true);
        add_hard_factor(new H1qqCorrectionA(), true);
        add_hard_factor(new H1qqCorrectionB(), true);
        add_hard_factor(new H1qqCorrection(), true);
        add_hard_factor(new H02gg(), true);
//         add_hard_factor(new H12gg(), true);
        add_hard_factor(new H12qqbar(), true);
        add_hard_factor(new H16gg(), true);
        add_hard_factor(new H1ggCorrection(), true);
//         add_hard_factor(new H112gq(), true);
//         add_hard_factor(new H122gq(), true);
        add_hard_factor(new H14gq(), true);
//         add_hard_factor(new H112qg(), true);
//         add_hard_factor(new H122qg(), true);
        add_hard_factor(new H14qg(), true);
    }
};

}

#endif // _HARDFACTORS_MOMENTUM_H_
