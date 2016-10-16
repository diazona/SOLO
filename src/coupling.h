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

// declared in integration/integrationcontext.h
class IntegrationContext;

typedef enum {KT, PT, Q1, Q2, Q3, KQ1, KQ2, KQ3} CouplingScale;

class Coupling {
public:
    Coupling() {}
    virtual ~Coupling() {}
    virtual double alphas(const IntegrationContext& ictx) const = 0;
    virtual const char* name() const = 0;
};

class FixedCoupling : public Coupling {
private:
    const double value;
    std::string _name;
public:
    FixedCoupling(const double alphas);
    double alphas(const IntegrationContext& ictx) const;
    const char* name() const;
};

class LORunningCoupling : public Coupling {
private:
    const double log_LambdaQCD;
    const double coefficient;
    const double regulator; // position of the Landau pole
    const CouplingScale m_scale_scheme;
    std::string _name;
public:
    LORunningCoupling(const double LambdaQCD, const double Ncbeta, const double regulator, const CouplingScale scale_scheme);
    double alphas(const IntegrationContext& ictx) const;
    const char* name() const;
};

std::ostream& operator<<(std::ostream& out, const Coupling& cpl);

#endif // _COUPLING_H_
