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

#ifndef _GLUONDIST_H_
#define _GLUONDIST_H_

#include <string>
#include "interp2d.h"

/**
 * A gluon distribution.
 */
class GluonDistribution {
public:
    /**
     * Return the value of the dipole gluon distribution at the given
     * values of r^2 and Q_s^2. r2 is the squared magnitude of the
     * dipole displacement vector x - y (equivalently, r2 is the
     * squared dipole size)
     */
    virtual double S2(double r2, double Qs2) = 0;
    /**
     * Return the value of the quadrupole gluon distribution at the given
     * values of r^2, s^2, t^2, and Q_s^2. r2 is the squared magnitude of
     * x - y, s2 is the squared magnitude of x - b, and t2 is the squared
     * magnitude of y - b, so that as vectors, r = s - t. r2, s2, and t2
     * thus are the squared lengths of three sides of a triangle.
     */
    virtual double S4(double r2, double s2, double t2, double Qs2) = 0;
    /**
     * Return the value of the momentum-space dipole gluon distribution
     * at the given values of q^2 and Q_s^2.
     */
    virtual double F(double q2, double Qs2) = 0;
    /**
     * Returns a human-readable name for the gluon distribution.
     */
    virtual const char* name() = 0;
};

/**
 * The GBW gluon distribution.
 */
class GBWGluonDistribution: public GluonDistribution {
public:
    /**
     * Returns the value of the GBW dipole gluon distribution,
     * exp(-r2 Qs2 / 4)
     */
    double S2(double r2, double Qs2);
    /**
     * Returns the value of the GBW quadrupole gluon distribution,
     * exp(-s2 Qs2 / 4) exp(-t2 Qs2 / 4)
     * This is just a product of two dipole distributions (valid
     * in the large-Nc limit)
     */
    double S4(double r2, double s2, double t2, double Qs2);
    /**
     * Returns the value of the GBW momentum space dipole gluon
     * distribution, exp(-q2 / Qs2) / (pi * Qs2)
     */
    double F(double q2, double Qs2);
    const char* name();
};

/**
 * The MV gluon distribution (aka "modified").
 * 
 * The position space dipole distribution takes the form
 * exp(-r2 Qs2 ln(e + 1 / (LambdaMV r)) / 4)
 * and the position space quadrupole distribution is computed
 * as a product of two factors of that form.
 * 
 * In momentum space, there is no analytic expression, so the
 * value of the distribution has to be determined numerically.
 * This implementation establishes a 2D grid in ln(q2) and ln(Qs2)
 * and computes the value of the distribution at the grid points
 * via a numerical integral. When F(q2, Qs2) is called, if the
 * values of q2 and Qs2 are within the boundaries of the grid,
 * the value of the distribution F is computed by 2D interpolation.
 * If q2 is smaller than the lower boundary of the grid, then
 * the value of the distribution F is computed from a series
 * expansion around q2 = 0. The series coefficients are
 * interpolated in Qs2 only.
 */
class MVGluonDistribution: public GluonDistribution {
public:
    /**
     * Constructs a new MV gluon distribution object.
     * 
     * `q2min`, `q2max`, `Qs2min`, and `Qs2max` specify the boundaries
     * of the region in which to interpolate the momentum space
     * distribution. The grid will be set up using these boundaries,
     * with a spacing automatically chosen to be reasonably accurate.
     * The series expansion used for `q2` < `q2min` is accurate to three
     * digits up to around `q2` = 1e-3 or so, but it's probably safer
     * to pass `q2min` around 1e-6. The other limits should be chosen
     * to include the range of values that will be needed.
     * 
     * `subinterval_limit` specifies the maximum number of subdivisions
     * used in computing the numeric integrals. Pass a larger value for
     * extreme parameters if the program crashes with a subdivision
     * error.
     */
    MVGluonDistribution(double LambdaMV, double q2min, double q2max, double Qs2min, double Qs2max, size_t subinterval_limit = 10000);
    ~MVGluonDistribution();
    
    /**
     * Returns the value of the MV dipole gluon distribution,
     * exp(-r2 Qs2 ln(e + 1 / (LambdaMV r)) / 4)
     */
    double S2(double r2, double Qs2);
    /**
     * Returns the value of the MV quadrupole gluon distribution,
     * computed as the product of two dipole gluon distributions,
     * S2(s2, Qs2) * S2 (t2, Qs2). This is valid in the large-Nc
     * limit.
     */
    double S4(double r2, double s2, double t2, double Qs2);
    /**
     * Returns the value of the MV momentum space dipole gluon
     * distribution.
     */
    double F(double q2, double Qs2);
    /**
     * Returns the name of the distribution, which incorporates
     * the values of the parameters.
     */
    const char* name();
private:
    double LambdaMV;
    double q2min, q2max;
    double Qs2min, Qs2max;
    /** Values of ln(q2) for the interpolation. */
    double* log_q2_values;
    /** Values of ln(Qs2) for the interpolation. */
    double* log_Qs2_values;
    /** Values of the leading coefficient in the series for small q2. */
    double* F_dist_leading_q2;
    /** Values of the subleading coefficient in the series for small q2. */
    double* F_dist_subleading_q2;
    /** Values of the gluon distribution for interpolation when q2 > q2min. */
    double* F_dist;
    
    // Structures used to do the interpolation
    gsl_interp* interp_dist_leading_q2;
    gsl_interp* interp_dist_subleading_q2;
    interp2d* interp_dist;
    gsl_interp_accel* q2_accel;
    gsl_interp_accel* Qs2_accel;
    
    std::string _name;
#ifdef GLUON_DIST_DRIVER
public:
    /**
     * This writes out the entire grid of the 2D interpolation to standard output.
     * 
     * Only used for testing.
     */
    void write_grid();
#endif
};

/** Prints the name of the gluon distribution to the given output stream. */
std::ostream& operator<<(std::ostream& out, GluonDistribution& gdist);

#endif // _GLUONDIST_H_