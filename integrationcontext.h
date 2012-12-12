#ifndef _INTEGRATIONCONTEXT_H_
#define _INTEGRATIONCONTEXT_H_

#include "context.h"

class IntegrationContext {
public:
    Context* ctx;
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
    double alphasbar;
    double qqfactor;
    double ggfactor;
    double gqfactor;
    double qgfactor;
    double S2r, S4rst;
    
    IntegrationContext(Context* ctx) :
      ctx(ctx),
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
      xp(0), xg(0),
      rx(0), ry(0), r2(0),
      sx(0), sy(0), s2(0),
      tx(0), ty(0), t2(0),
      Qs2(0),
      alphasbar(0),
      qqfactor(0),
      ggfactor(0),
      gqfactor(0),
      qgfactor(0),
      S2r(0), S4rst(0) {
    };
    void update_positions(double xx, double xy, double yx, double yy, double bx, double by);
    void update_momenta(double q1x, double q1y, double q2x, double q2y, double q3x, double q3y);
    void update_auxiliary(double xiprime);
    void update_parton_factors(double z, double y);
    void set_xi_to_1() { update_parton_factors(z, 1); }
};

#endif // _INTEGRATIONCONTEXT_H_
