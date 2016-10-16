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

#include <sstream>
#include "integration/integrationcontext.h"
#include "coupling.h"

using namespace std;


static inline double square(double x, double y) {
    return x * x + y * y;
}

static inline double scale(const IntegrationContext& ictx, const CouplingScale scale) {
    switch (scale) {
        case KT:
            return ictx.kT2;
        case PT:
            return ictx.ctx.pT2;
        case Q1:
            return ictx.q12;
        case Q2:
            return ictx.q22;
        case Q3:
            return ictx.q32;
        case KQ1:
            return square(ictx.kT - ictx.q1x, ictx.q1y);
        case KQ2:
            return square(ictx.kT - ictx.q2x, ictx.q2y);
        case KQ3:
            return square(ictx.kT - ictx.q3x, ictx.q3y);
        default:
            assert(false);
    }
}


FixedCoupling::FixedCoupling(const double alphas) : value(alphas) {
    ostringstream s;
    s << "Fixed(alphas = " << value << ")";
    _name = s.str();
}
double FixedCoupling::alphas(const IntegrationContext& ictx) const {
    return value;
}
const char* FixedCoupling::name() const {
    return _name.c_str();
}

LORunningCoupling::LORunningCoupling(const double LambdaQCD, const double Ncbeta, const double regulator, const CouplingScale scale_scheme) :
  log_LambdaQCD(log(LambdaQCD)),
  coefficient(M_PI / Ncbeta), // π/(Nc × β)
  regulator(regulator),
  m_scale_scheme(scale_scheme) {
    ostringstream s;
    s << "LORunning(LambdaQCD = " << LambdaQCD << ", Nc × β = " << Ncbeta << ", regulator = " << regulator << ")";
    _name = s.str();
}
double LORunningCoupling::alphas(const IntegrationContext& ictx) const {
    double l_scale = scale(ictx, m_scale_scheme);
    return coefficient / (log(l_scale + regulator) - log_LambdaQCD);
}
const char* LORunningCoupling::name() const {
    return _name.c_str();
}


ostream& operator<<(ostream& out, const Coupling& cpl) {
    out << cpl.name();
    return out;
}
