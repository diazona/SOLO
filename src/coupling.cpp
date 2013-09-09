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
#include "coupling.h"

using namespace std;

FixedCoupling::FixedCoupling(double alphas) : value(alphas) {
    ostringstream s;
    s << "Fixed(alphas = " << value << ")";
    _name = s.str();
}
double FixedCoupling::alphas(double kT2) {
    return value;
}
const char* FixedCoupling::name() {
    return _name.c_str();
}

LORunningCoupling::LORunningCoupling(double LambdaQCD, double Ncbeta, double regulator) :
    log_LambdaQCD(log(LambdaQCD)),
    coefficient(M_PI / Ncbeta), // π/(Nc × β)
    regulator(regulator) {
    ostringstream s;
    s << "LORunning(LambdaQCD = " << LambdaQCD << ", Nc × β = " << Ncbeta << ", regulator = " << regulator << ")";
    _name = s.str();
}
double LORunningCoupling::alphas(double kT2) {
    return coefficient / (log(kT2 + regulator) - log_LambdaQCD);
}
const char* LORunningCoupling::name() {
    return _name.c_str();
}


std::ostream& operator<<(std::ostream& out, Coupling& cpl) {
    out << cpl.name();
    return out;
}
