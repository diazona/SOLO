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

#ifndef _COUPLING_H_
#define _COUPLING_H_

#include <cmath>
#include <ostream>
#include <string>

class Coupling {
public:
    Coupling() {}
    virtual ~Coupling() {}
    virtual double alphas(double kT2) = 0;
    virtual const char* name() = 0;
};

class FixedCoupling : public Coupling {
private:
    double value;
    std::string _name;
public:
    FixedCoupling(double alphas);
    double alphas(double kT2);
    const char* name();
};

class LORunningCoupling : public Coupling {
private:
    double log_LambdaQCD;
    double coefficient;
    double regulator; // position of the Landau pole
    std::string _name;
public:
    LORunningCoupling(double LambdaQCD, double Ncbeta, double regulator);
    double alphas(double kT2);
    const char* name();
};

std::ostream& operator<<(std::ostream& out, Coupling& cpl);

#endif // _COUPLING_H_
