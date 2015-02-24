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

#include "hardfactor.h"

namespace position {
    const PositionIntegrationType dipole(2);
    const PositionIntegrationType quadrupole(4);
}

namespace radial {
    const AngleIndependentPositionIntegrationType dipole(1);
    const AngleIndependentPositionIntegrationType quadrupole(2);
    const RescaledAngleIndependentPositionIntegrationType rescaled_dipole(1);
    const RescaledAngleIndependentPositionIntegrationType rescaled_quadrupole(2);
}

namespace momentum {
    const NoIntegrationType none;
    const MomentumIntegrationType momentum1(2);
    const MomentumIntegrationType momentum2(4);
    const MomentumIntegrationType momentum3(6);
    const RadialMomentumIntegrationType radialmomentum1(2);
    const RadialMomentumIntegrationType radialmomentum2(4);
    const RadialMomentumIntegrationType radialmomentum3(6);
    const XiPIntegrationType momentumxip1(3);
    const XiPIntegrationType momentumxip2(5);
    const QLimitedMomentumIntegrationType qlim(2);
}

// for weird technical reasons, this needs to be in a .cpp file, not a header
HardFactor::~HardFactor() {}
