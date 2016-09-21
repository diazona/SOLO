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

#ifndef _INTEGRATIONCONTEXT_H_
#define _INTEGRATIONCONTEXT_H_

#include "../configuration/context.h"

class Modifiers {
public:
    bool exact_xg;
    bool divide_xi;

    Modifiers() : exact_xg(false), divide_xi(true) {}
    bool operator<(const Modifiers& other) const;
    bool operator==(const Modifiers& other) const;
};

class IntegrationContext {
public:
    const Context& ctx;
    const ThreadLocalContext& tlctx;
    // updated
    double z;
    double xi;
    double xx, xy;
    double yx, yy;
    double bx, by;
    // calculated
    double z2, xi2;
    double kT2, kT;
    double xp, xg;
    double Yg;
    double rx, ry, r2;
    double sx, sy, s2;
    double tx, ty, t2;
    // some of these are updated
    double q1x, q1y, q12;
    double q2x, q2y, q22;
    double q3x, q3y, q32;
    double Qs2;
    double mu2;
    double alphas, alphas_2pi; // αs/2π
    double qqfactor;
    double ggfactor;
    double gqfactor;
    double qgfactor;
    double S2r, S4rst;
    double Fk;
    double Fq1, Fq2, Fq3;
    double Fkq1, Fkq2, Fkq3;

    IntegrationContext(const Context& ctx, const ThreadLocalContext& tlctx) :
      ctx(ctx),
      tlctx(tlctx),
      z(0), xi(0),
      xx(0), xy(0),
      yx(0), yy(0),
      bx(0), by(0),
      q1x(0), q1y(0),
      q2x(0), q2y(0),
      q3x(0), q3y(0),
      z2(0), xi2(0),
      kT2(0), kT(0),
      xp(0), xg(0),
      rx(0), ry(0), r2(0),
      sx(0), sy(0), s2(0),
      tx(0), ty(0), t2(0),
      Qs2(0),
      mu2(0),
      alphas(0),
      alphas_2pi(0),
      qqfactor(0),
      ggfactor(0),
      gqfactor(0),
      qgfactor(0),
      S2r(0), S4rst(0),
      Fk(0),
      Fq1(0), Fq2(0), Fq3(0),
      Fkq1(0), Fkq2(0), Fkq3(0) {
    };

    void recalculate_everything(const bool exact_xg, const bool divide_xi);
    void recalculate_everything_from_position(const bool quadrupole, const bool divide_xi);
    void recalculate_everything_from_momentum(const size_t dimensions, const bool exact_xg, const bool divide_xi);

private:
    void recalculate_from_position(const bool quadrupole);
    void recalculate_from_momentum(const size_t dimensions);
    void recalculate_longitudinal(const bool exact_xg);
    void recalculate_coupling();
    void recalculate_position_gdist(const bool quadrupole);
    void recalculate_momentum_gdist(const size_t dimensions);
    void recalculate_parton_functions(const bool divide_xi);
};

#endif // _INTEGRATIONCONTEXT_H_
