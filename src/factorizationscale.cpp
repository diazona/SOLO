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

#include <cassert>
#include <sstream>
#include <gsl/gsl_math.h>
#include "factorizationscale.h"
#include "integration/integrationcontext.h"

using namespace std;

FixedFactorizationScale::FixedFactorizationScale(double mu2) : value(mu2) {
    ostringstream s;
    s << "Fixed(mu2 = " << value << ")";
    _name = s.str();
}
double FixedFactorizationScale::mu2(const IntegrationContext& ictx) {
    return value;
}
const char* FixedFactorizationScale::name() {
    return _name.c_str();
}

PTProportionalFactorizationScale::PTProportionalFactorizationScale(double coefficient) : coefficient(coefficient) {
    ostringstream s;
    s << "PTProportional(coefficient = " << coefficient << ")";
    _name = s.str();
}
double PTProportionalFactorizationScale::mu2(const IntegrationContext& ictx) {
    return coefficient * ictx.ctx.pT2;
}
const char* PTProportionalFactorizationScale::name() {
    return _name.c_str();
}

RPerpFactorizationScale::RPerpFactorizationScale(double coefficient) : coefficient(coefficient) {
    ostringstream s;
    s << "RPerp(coefficient = " << coefficient << ")";
    _name = s.str();
}
double RPerpFactorizationScale::mu2(const IntegrationContext& ictx) {
    assert(ictx.r2 > 0);
    // implement the CSS r regularization
    if (ictx.ctx.css_r_regularization) {
        return coefficient * (1 + ictx.r2 / ictx.ctx.css_r2_max) / ictx.r2;
    }
    else {
        return coefficient / ictx.r2;
    }
}
const char* RPerpFactorizationScale::name() {
    return _name.c_str();
}



std::ostream& operator<<(std::ostream& out, FactorizationScale& fs) {
    out << fs.name();
    return out;
}
