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

#ifndef _HARDFACTORS_POSITION_H_
#define _HARDFACTORS_POSITION_H_

// Hard factors in position space

#include "hardfactor.h"
#include "../integration/integrationtype.h"

namespace position {

extern const PositionIntegrationType dipole;
extern const PositionIntegrationType quadrupole;

class H02qq : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H02qq";
    }
    const IntegrationType* get_type() const {
        return &dipole;
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H12qq : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H12qq";
    }
    const IntegrationType* get_type() const {
        return &dipole;
    }
    void Fs(const IntegrationContext* ictx, double* real, double* imag) const;
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H14qqPrimary : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H14qq";
    }
    const IntegrationType* get_type() const {
        return &quadrupole;
    }
    void Fs(const IntegrationContext* ictx, double* real, double* imag) const;
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H14qqResidual : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H14qqResidual";
    }
    const IntegrationType* get_type() const {
        return &dipole;
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H14qq : public HardFactor {
public:
    H14qq() {
        terms[0] = new H14qqPrimary();
        terms[1] = new H14qqResidual();
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

class H02gg : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H02gg";
    }
    const IntegrationType* get_type() const {
        return &dipole;
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H12gg : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H12gg";
    }
    const IntegrationType* get_type() const {
        return &dipole;
    }
    void Fs(const IntegrationContext* ictx, double* real, double* imag) const;
    void Fn(const IntegrationContext* ictx, double* real, double* imag) const;
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H12qqbarPrimary : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H12qqbarPrimary";
    }
    const IntegrationType* get_type() const {
        return &quadrupole;
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H12qqbarResidual : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H12qqbarResidual";
    }
    const IntegrationType* get_type() const {
        return &dipole;
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H12qqbar : public HardFactor {
public:
    H12qqbar() {
        terms[0] = new H12qqbarPrimary();
        terms[1] = new H12qqbarResidual();
    }
    ~H12qqbar() {
        delete terms[0];
        terms[0] = NULL;
        delete terms[1];
        terms[1] = NULL;
    }
    const char* get_name() const {
        return "H12qqbar";
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

class H16ggPrimary : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H16ggPrimary";
    }
    const IntegrationType* get_type() const {
        return &quadrupole;
    }
    void Fs(const IntegrationContext* ictx, double* real, double* imag) const;
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H16ggResidual : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H16ggResidual";
    }
    const IntegrationType* get_type() const {
        return &dipole;
    }
    void Fd(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H16gg : public HardFactor {
public:
    H16gg() {
        terms[0] = new H16ggPrimary();
        terms[1] = new H16ggResidual();
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

class H112gq : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H112gq";
    }
    const IntegrationType* get_type() const {
        return &dipole;
    }
    void Fn(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H122gq : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H122gq";
    }
    const IntegrationType* get_type() const {
        return &dipole;
    }
    void Fn(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H14gq : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H14gq";
    }
    const IntegrationType* get_type() const {
        return &quadrupole;
    }
    void Fn(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H112qg : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H112qg";
    }
    const IntegrationType* get_type() const {
        return &dipole;
    }
    void Fn(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H122qg : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H122qg";
    }
    const IntegrationType* get_type() const {
        return &dipole;
    }
    void Fn(const IntegrationContext* ictx, double* real, double* imag) const;
};

class H14qg : public HardFactorTerm {
public:
    const char* get_name() const {
        return "H14qg";
    }
    const IntegrationType* get_type() const {
        return &quadrupole;
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
    registry() {
        add_hard_factor(new H02qq(), true);
        add_hard_factor(new H12qq(), true);
        add_hard_factor(new H14qq(), true);
        add_hard_factor(new H02gg(), true);
        add_hard_factor(new H12gg(), true);
        add_hard_factor(new H12qqbar(), true);
        add_hard_factor(new H16gg(), true);
        add_hard_factor(new H112gq(), true);
        add_hard_factor(new H122gq(), true);
        add_hard_factor(new H14gq(), true);
        add_hard_factor(new H112qg(), true);
        add_hard_factor(new H122qg(), true);
        add_hard_factor(new H14qg(), true);
    }
};

}

#endif // _HARDFACTORS_POSITION_H_