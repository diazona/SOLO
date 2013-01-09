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

FixedCoupling::FixedCoupling(double alphasbar) : value(alphasbar) {
    ostringstream s;
    s << "Fixed(alphasbar = " << value << ")";
    _name = s.str();
}
double FixedCoupling::alphasbar(double kT2) {
    return value;
}
const char* FixedCoupling::name() {
    return _name.c_str();
}

LORunningCoupling::LORunningCoupling(double LambdaQCD, double beta, double regulator) :
    log_LambdaQCD(log(LambdaQCD)),
    inverse_beta_2pi(0.5 / (M_PI * beta)),
    regulator(regulator) {
    ostringstream s;
    s << "LORunning(LambdaQCD = " << LambdaQCD << ", beta = " << beta << ", regulator = " << regulator << ")";
    _name = s.str();
}
double LORunningCoupling::alphasbar(double kT2) {
    return inverse_beta_2pi / (log(kT2 + regulator) - log_LambdaQCD);
}
const char* LORunningCoupling::name() {
    return _name.c_str();
}


std::ostream& operator<<(std::ostream& out, Coupling& cpl) {
    out << cpl.name();
    return out;
}
