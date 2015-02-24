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

#ifndef _HARDFACTORS_RADIAL_H_
#define _HARDFACTORS_RADIAL_H_

// Hard factors in position space with radial integration

#include "hardfactor.h"
#include "../integration/integrationtype.h"

namespace radial {

extern const AngleIndependentPositionIntegrationType dipole;
extern const AngleIndependentPositionIntegrationType quadrupole;
extern const RescaledAngleIndependentPositionIntegrationType rescaled_dipole;
extern const RescaledAngleIndependentPositionIntegrationType rescaled_quadrupole;

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

#endif // _HARDFACTORS_RADIAL_H_
