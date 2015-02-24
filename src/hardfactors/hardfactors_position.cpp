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

// Hard factors in position space

#include <cassert>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include "hardfactor.h"
#include "hardfactors_position.h"

#define checkfinite(d) assert(gsl_finite(d))

namespace position {
    const PositionIntegrationType dipole(2);
    const PositionIntegrationType quadrupole(4);
}
