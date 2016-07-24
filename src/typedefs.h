/*
 * Part of SOLO
 *
 * Copyright 2014 David Zaslavsky
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

#pragma once

#include <map>
#include <vector>
#include "integration/integrationregion.h"

// declared in hardfactor.h
class HardFactor;
class HardFactorTerm;

typedef std::vector<const HardFactor*> HardFactorList;
typedef std::vector<const HardFactorTerm*> HardFactorTermList;
typedef std::map<const IntegrationRegion*, HardFactorTermList, bool(*)(const IntegrationRegion*, const IntegrationRegion*)> HardFactorTypeMap;
