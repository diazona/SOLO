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
    double xp, xg, xa;
    double Yg, Ya;
    double rx, ry, r2;
    double sx, sy, s2;
    double tx, ty, t2;
    // some of these are updated
    double q1x, q1y, q12;
    double q2x, q2y, q22;
    double q3x, q3y, q32;
    double xiprime;
    double xiprime2;
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
    double qmax;

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
      xiprime(0),
      xiprime2(0),
      z2(0), xi2(0),
      kT2(0), kT(0),
      xp(0), xg(0), xa(0),
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
      Fkq1(0), Fkq2(0), Fkq3(0),
      qmax(0) {
    };
    void update_kinematics(double z, double xi, size_t core_dimensions);
    void update_positions(double xx, double xy, double yx, double yy, double bx, double by);
    void update_momenta(double q1x, double q1y, double q2x, double q2y, double q3x, double q3y);
    void update_auxiliary(double xiprime);
    void update_parton_functions();
    void set_xi_to_1(size_t core_dimensions) {
        update_kinematics(z, 1, core_dimensions);
        update_parton_functions();
    }
};

#endif // _INTEGRATIONCONTEXT_H_
