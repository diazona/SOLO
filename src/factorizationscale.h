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

#pragma once

#ifndef _FACTORIZATIONSCALE_H_
#define _FACTORIZATIONSCALE_H_

#include <cmath>
#include <ostream>
#include <string>

// declared in integrationcontext.h
class IntegrationContext;

class FactorizationScale {
public:
    virtual double mu2(const IntegrationContext& ictx) = 0;
    virtual const char* name() = 0;
};

class FixedFactorizationScale : public FactorizationScale {
private:
    double value;
    std::string _name;
public:
    FixedFactorizationScale(double mu2);
    double mu2(const IntegrationContext& ictx);
    const char* name();
};

class PTProportionalFactorizationScale : public FactorizationScale {
private:
    double coefficient;
    std::string _name;
public:
    PTProportionalFactorizationScale(double coefficient);
    double mu2(const IntegrationContext& ictx);
    const char* name();
};

class RPerpFactorizationScale : public FactorizationScale {
private:
    double coefficient;
    std::string _name;
public:
    RPerpFactorizationScale(double coefficient);
    double mu2(const IntegrationContext& ictx);
    const char* name();
};

std::ostream& operator<<(std::ostream& out, FactorizationScale& fs);

#endif // _FACTORIZATIONSCALE_H_
