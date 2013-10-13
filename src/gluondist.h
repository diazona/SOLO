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

#include <exception>
#include <string>
#include "interp2d.h"

class NoGridException : public std::exception {
};

class GluonDistributionS2RangeException : public std::exception {
private:
    double _r2, _Y;
protected:
    std::string _message;
public:
    GluonDistributionS2RangeException(double r2, double Y) throw() :
        _r2(r2), _Y(Y) {
        std::ostringstream s;
        s << "Gluon distribution S2 evaluated at " << _r2 << "," << Y << " (out of range)";
        _message = s.str();
    }
    GluonDistributionS2RangeException(const GluonDistributionS2RangeException& e) throw() :
        _r2(e._r2), _Y(e._Y), _message(e._message) {
    }
    GluonDistributionS2RangeException& operator=(const GluonDistributionS2RangeException& e) throw() {
        _r2 = e._r2;
        _Y = e._Y;
        _message = e._message;
    }
    virtual ~GluonDistributionS2RangeException() throw() {
    }
    const double r2() const {
        return _r2;
    }
    const double Y() const {
        return _Y;
    }
    const char* what() const throw() {
        return _message.c_str();
    }
};

class GluonDistributionS4RangeException : public GluonDistributionS2RangeException {
private:
    double _s2, _t2;
public:
    GluonDistributionS4RangeException(double r2, double s2, double t2, double Y) throw() :
        GluonDistributionS2RangeException(r2, Y), _s2(s2), _t2(t2) {
        std::ostringstream s;
        s << "Gluon distribution S4 evaluated at " << this->r2() << "," << this->s2() << "," << this->t2() << "," << this->Y() << " (out of range)";
        _message = s.str();
    }
    GluonDistributionS4RangeException(const GluonDistributionS4RangeException& e) throw() :
        GluonDistributionS2RangeException(e), _s2(e._s2), _t2(e._t2) {
    }
    GluonDistributionS4RangeException& operator=(const GluonDistributionS4RangeException& e) throw() {
        GluonDistributionS2RangeException::operator=(e);
        _s2 = e._s2;
        _t2 = e._t2;
    }
    virtual ~GluonDistributionS4RangeException() throw() {
    }
    const double s2() const {
        return _s2;
    }
    const double t2() const {
        return _t2;
    }
};

class GluonDistributionFRangeException : public std::exception {
private:
    double _q2, _Y;
    std::string _message;
public:
    GluonDistributionFRangeException(double q2, double Y) throw() :
        _q2(q2), _Y(Y) {
        std::ostringstream s;
        s << "Gluon distribution F evaluated at " << _q2 << "," << _Y << " (out of range)";
        _message = s.str();
    }
    GluonDistributionFRangeException(const GluonDistributionFRangeException& e) throw() :
        _q2(e._q2), _Y(e._Y), _message(e._message) {
    }
    GluonDistributionFRangeException& operator=(const GluonDistributionFRangeException& e) throw() {
        _q2 = e._q2;
        _Y = e._Y;
        _message = e._message;
    }
    virtual ~GluonDistributionFRangeException() throw() {
    }
    const double q2() const {
        return _q2;
    }
    const double Y() const {
        return _Y;
    }
    const char* what() const throw() {
        return _message.c_str();
    }
};

/**
 * A gluon distribution.
 */
class GluonDistribution {
public:
    GluonDistribution() {};
    virtual ~GluonDistribution() {};
    /**
     * Return the value of the dipole gluon distribution at the given
     * values of r^2 and Y. r2 is the squared magnitude of the
     * dipole displacement vector x - y (equivalently, r2 is the
     * squared dipole size).
     */
    virtual double S2(double r2, double Y) = 0;
    /**
     * Return the value of the quadrupole gluon distribution at the given
     * values of r^2, s^2, t^2, and Y. r2 is the squared magnitude of
     * x - y, s2 is the squared magnitude of x - b, and t2 is the squared
     * magnitude of y - b, so that as vectors, r = s - t. r2, s2, and t2
     * thus are the squared lengths of three sides of a triangle.
     */
    virtual double S4(double r2, double s2, double t2, double Y) = 0;
    /**
     * Return the value of the momentum-space dipole gluon distribution
     * at the given values of q^2 and Y.
     */
    virtual double F(double q2, double Y) = 0;
    /**
     * Return the saturation scale corresponding to the given value of Y.
     */
    virtual double Qs2(const double Y) const = 0;
    /**
     * Returns a human-readable name for the gluon distribution.
     */
    virtual const char* name() = 0;
    
    /**
     * This writes out the entire grid of the 2D momentum interpolation to standard output.
     * 
     * Only used for testing.
     */
    virtual void write_pspace_grid(std::ostream& out);
    /**
     * This writes out the entire grid of the 2D position interpolation to standard output.
     * 
     * Only used for testing.
     */
    virtual void write_rspace_grid(std::ostream& out);
    /**
     * This writes out the entire grid of the 1D saturation scale interpolation to standard output.
     * 
     * Only used for testing.
     */
    virtual void write_satscale_grid(std::ostream& out);
};

/**
 * The GBW gluon distribution.
 */
class GBWGluonDistribution: public GluonDistribution {
public:
    GBWGluonDistribution(double Q02, double x0, double lambda);
    /**
     * Returns the value of the GBW dipole gluon distribution,
     * exp(-r2 Qs2 / 4)
     */
    double S2(double r2, double Y);
    /**
     * Returns the value of the GBW quadrupole gluon distribution,
     * exp(-s2 Qs2 / 4) exp(-t2 Qs2 / 4)
     * This is just a product of two dipole distributions (valid
     * in the large-Nc limit)
     */
    double S4(double r2, double s2, double t2, double Y);
    /**
     * Returns the value of the GBW momentum space dipole gluon
     * distribution, exp(-q2 / Qs2) / (pi * Qs2)
     */
    double F(double q2, double Y);
    /**
     * Returns the standard saturation scale, Q0^2(x0/x)^λ
     */
    double Qs2(const double Y) const;
    const char* name();
private:
    double Q02x0lambda;
    double lambda;
    std::string _name;
};

/**
 * An abstract base class for gluon distributions which only have a position
 * space or momentum space definition - that is, for when you have an analytic
 * expression for S(r2, Qs2) or F(q2, Qs2), but not both. The values in the
 * other space are computed automatically via a numerical integral.
 * 
 * This is not meant to be subclassed directly. Gluon distributions should be
 * subclasses of either AbstractPositionGluonDistribution or
 * AbstractMomentumGluonDistribution. This just encapsulates the behavior
 * common to both.
 * 
 * This implementation establishes a 2D grid in ln(u2) and ln(Qs2), where u
 * represents either r or q, and numerically integrates to get the values of
 * the gluon distribution at the grid points. When F(q2, Qs2) or S2(r2, Qs2) is
 * called, if the values of u2 and Qs2 are within the boundaries of the grid,
 * the value of the distribution S2 or F is computed by 2D interpolation. If u2
 * is smaller than the lower boundary of the grid, then the value of the gluon
 * distribution is computed from a series expansion around u2 = 0. The series
 * coefficients are interpolated in Qs2 only.
 */
class AbstractTransformGluonDistribution : public GluonDistribution {
public:
    /**
     * Constructs a new transform gluon distribution object.
     * 
     * `u2min`, `u2max`, `Ymin`, and `Ymax` specify the boundaries
     * of the region in which to interpolate the momentum space
     * distribution. The grid will be set up using these boundaries,
     * with a spacing automatically chosen to be reasonably accurate.
     * The series expansion used for `u2` < `u2min` is accurate to three
     * digits up to around `u2` = 1e-3 or so, but it's probably safer
     * to pass `u2min` around 1e-6. The other limits should be chosen
     * to include the range of values that will be needed.
     * 
     * `subinterval_limit` specifies the maximum number of subdivisions
     * used in computing the numeric integrals. Pass a larger value for
     * extreme parameters if the program crashes with a subdivision
     * error.
     */
    AbstractTransformGluonDistribution(double u2min, double u2max, double Ymin, double Ymax, size_t subinterval_limit, double (*gdist_integrand)(double, void*), double (*gdist_series_term_integrand)(double, void*));
    virtual ~AbstractTransformGluonDistribution();

    /**
     * Returns the value of the position space quadrupole gluon distribution,
     * by default computed as the product of two dipole gluon distributions,
     * S2(s2, Qs2) * S2 (t2, Qs2). This is usually valid in the large-Nc
     * limit.
     */
    virtual double S4(double r2, double s2, double t2, double Y);

    /**
     * Returns the name of the distribution, which should incorporate
     * the values of the parameters.
     */
    virtual const char* name() = 0;

protected:
    /**
     * Handles the actual calculation of the points to use for interpolation.
     * This should be called from the constructor of each subclass that
     * implements F or S2.
     */
    void setup();
    
    /**
     * Returns the value of the dipole distribution, F or S2.
     */
    double dipole_distribution(double u2, double Y);
    /**
     * Write a representation of the internal interpolation grid to the given
     * output stream.
     */
    void write_grid(std::ostream& out);
    
private:
    double u2min, u2max;
    double Ymin, Ymax;
    /** Values of ln(u2) for the interpolation. */
    double* log_u2_values;
    /** Values of Y for the interpolation. */
    double* Y_values;
    
    // G stands for either F or S2
    /** Values of the leading coefficient in the series for small u2. */
    double* G_dist_leading_u2;
    /** Values of the subleading coefficient in the series for small u2. */
    double* G_dist_subleading_u2;
    /** Values of the gluon distribution for interpolation when u2 > u2min. */
    double* G_dist;

    // The functions used to compute the transformed distribution
    double (*gdist_integrand)(double, void*);
    double (*gdist_series_term_integrand)(double, void*);
    
    // Structures used to do the interpolation
    gsl_interp* interp_dist_leading_u2;
    gsl_interp* interp_dist_subleading_u2;
    // one of interp_dist_1D or interp_dist_2D will be NULL and the other one
    // will be used, depending on whether there is a range of Y values or
    // just a single one
    gsl_interp* interp_dist_1D;
    interp2d* interp_dist_2D;
    
    gsl_interp_accel* u2_accel;
    gsl_interp_accel* Y_accel;
    
    size_t u2_dimension;
    size_t Y_dimension;
    
    size_t subinterval_limit;
};

/**
 * An abstract base class for gluon distributions which only have an analytic
 * definition in position space. The momentum space distribution is computed
 * automatically via a numerical integral. See the documentation for
 * AbstractTransformGluonDistribution for details.
 */
class AbstractPositionGluonDistribution : public AbstractTransformGluonDistribution {
public:
    /**
     * Constructs a new position space gluon distribution object.
     * 
     * `q2min`, `q2max`, `Ymin`, and `Ymax` specify the boundaries
     * of the region in which to interpolate the momentum space
     * distribution in two dimensions. Calls to F(q2, Y) will succeed
     * for `q2 <= q2max` and `Ymin <= Y <= Ymax`.
     */
    AbstractPositionGluonDistribution(double q2min, double q2max, double Ymin, double Ymax, size_t subinterval_limit = 10000);
    virtual ~AbstractPositionGluonDistribution() {};
    
    /**
     * Returns the value of the momentum space dipole gluon distribution.
     * 
     * If `q2min <= q2 <= q2max` and `Ymin <= Y <= Ymax` (where `q2min`,
     * `q2max`, `Ymin`, and `Ymax` are the values passed to the constructor),
     * the result for `F` is computed using a 2D interpolation. Otherwise, if
     * `q2 < q2min` and `Ymin <= Y <= Ymax`, it is computed using a two-term
     * series expansion in `q2` with coefficients obtained by 1D interpolation
     * in `Y`. Otherwise, this throws a GluonDistributionFRangeException.
     */
    double F(double q2, double Y);
    /**
     * Write a representation of the internal momentum space grid to the
     * given output stream.
     */
    void write_pspace_grid(std::ostream& out);
};

/**
 * An abstract base class for gluon distributions which only have an analytic
 * definition in momentum space. The position space distribution is computed
 * automatically via a numerical integral. See the documentation for
 * AbstractTransformGluonDistribution for details.
 */
class AbstractMomentumGluonDistribution : public AbstractTransformGluonDistribution {
public:
    /**
     * Constructs a new momentum space gluon distribution object.
     * 
     * `r2min`, `r2max`, `Ymin`, and `Ymax` specify the boundaries
     * of the region in which to interpolate the momentum space
     * distribution in two dimensions. Calls to S2(r2, Y) will succeed
     * for `r2 <= r2max` and `Ymin <= Y <= Ymax`.
     */
    AbstractMomentumGluonDistribution(double r2min, double r2max, double Ymin, double Ymax, size_t subinterval_limit = 10000);
    virtual ~AbstractMomentumGluonDistribution() {};

    /**
     * Returns the value of the position space dipole gluon distribution.
     * 
     * If `r2min <= r2 <= r2max` and `Ymin <= Y <= Ymax` (where `r2min`,
     * `r2max`, `Ymin`, and `Ymax` are the values passed to the constructor),
     * the result for `S2` is computed using a 2D interpolation. Otherwise, if
     * `r2 < r2min` and `Ymin <= Y <= Ymax`, it is computed using a two-term
     * series expansion in `r2` with coefficients obtained by 1D interpolation
     * in `Y`. Otherwise, this throws a GluonDistributionS2RangeException.
     */
    double S2(double r2, double Y);

    /**
     * Write a representation of the internal position space grid to the
     * given output stream.
     */
    void write_rspace_grid(std::ostream& out);
};

/**
 * The MV gluon distribution.
 * 
 * The position space dipole distribution takes the form
 * exp(-(r2 Qs02MV)^gammaMV ln(e + 1 / (LambdaMV r)) / 4)
 * and the position space quadrupole distribution is computed
 * as a product of two factors of that form.
 */
class MVGluonDistribution: public AbstractPositionGluonDistribution {
public:
    /**
     * Constructs a new MV gluon distribution object.
     */
    MVGluonDistribution(
        double LambdaMV,
        double gammaMV,
        double q2min,
        double q2max,
        double Ymin,
        double Ymax,
        double Q02,
        double x0,
        double lambda,
        size_t subinterval_limit = 10000);
    virtual ~MVGluonDistribution() {};
    
    /**
     * Returns the value of the MV dipole gluon distribution,
     * exp(-(r2 Qs02MV)^gammaMV ln(e + 1 / (LambdaMV r)) / 4)
     */
    double S2(double r2, double Y);
    /**
     * Returns the standard saturation scale, Q0^2(x0/x)^λ
     */
    double Qs2(const double Y) const;
    /**
     * Returns the name of the distribution, which incorporates
     * the values of the parameters.
     */
    const char* name();
protected:
    double LambdaMV;
    double gammaMV;
    
    std::string _name;
private:
    double Q02x0lambda;
    double lambda;
};

/**
 * A modified version of the MV gluon distribution which ignores the actual
 * saturation scale and uses a fixed value instead.
 */
class FixedSaturationMVGluonDistribution : public MVGluonDistribution {
public:
    /**
     * Constructs a new modified gluon distribution object.
     */
    FixedSaturationMVGluonDistribution(double LambdaMV, double gammaMV, double q2min, double q2max, double YMV, double Q02, double x0, double lambda, size_t subinterval_limit = 10000);
    virtual ~FixedSaturationMVGluonDistribution() {};
    
    /**
     * Returns the value of the dipole gluon distribution,
     * exp(-(r2 Qs02)^gammaMV ln(e + 1 / (LambdaMV r)) / 4)
     * 
     * The parameter Y is not used.
     */
    double S2(double r2, double Y);
protected:
    double YMV;
};

/**
 * An unintegrated gluon distribution that is flat below Qs and has a power
 * tail above it.
 */
class PlateauPowerGluonDistribution : public AbstractMomentumGluonDistribution {
public:
    PlateauPowerGluonDistribution(
        double gamma,
        double r2min,
        double r2max,
        double Ymin,
        double Ymax,
        double Q02,
        double x0,
        double lambda,
        size_t subinterval_limit = 10000);
    virtual ~PlateauPowerGluonDistribution();
    /**
     * Returns the value of the dipole gluon distribution,
     * exp(-(r2 Qs02MV)^gammaMV ln(e + 1 / (LambdaMV r)) / 4)
     */
    double F(double q2, double Y);
    /**
     * Returns the standard saturation scale, Q0^2(x0/x)^λ
     */
    double Qs2(const double Y) const;
    /**
     * Returns the name of the distribution, which incorporates
     * the values of the parameters.
     */
    const char* name();
protected:
    double gammaPP;
    
    std::string _name;
private:
    double Q02x0lambda;
    double lambda;
};

class FileDataGluonDistribution : public GluonDistribution {
public:
    typedef enum satscale_source {NONE, POSITION_THRESHOLD, MOMENTUM_THRESHOLD} satscale_source_type;
    /**
     * Constructs a new gluon distribution reading from the specified file.
     */
    FileDataGluonDistribution(std::string pos_filename, std::string mom_filename, double xinit, satscale_source_type satscale_source, double satscale_threshold);
    /**
     * Constructs a new gluon distribution reading from the specified file.
     * This constructor uses the default saturation scale.
     */
    FileDataGluonDistribution(std::string pos_filename, std::string mom_filename, double Q02, double x0, double lambda, double xinit);
    virtual ~FileDataGluonDistribution();
    
    /**
     * Returns the interpolated value of the gluon distribution at the specified values.
     */
    double S2(double r2, double Y);
    
    double S4(double r2, double s2, double t2, double Y);
    
    double Qs2(const double Y) const;
    
    double F(double q2, double Y);
    
    /**
     * Returns the derivative of F with respect to q2 at the given
     * values of q2 and Y.
     */
    double Fprime(double q2, double Y);
    
    /**
     * Returns the second derivative of F with respect to q2 at the given
     * values of q2 and Y.
     */
    double Fpprime(double k2, double Y);
    
    /**
     * Returns the derivative of S with respect to r2 at the given
     * values of r2 and Y.
     */
    double S2prime(double r2, double Y);
    
    /**
     * Returns the second derivative of S2 with respect to r2 at the given
     * values of r2 and Y.
     */
    double S2pprime(double r2, double Y);
    
    /**
     * Returns the name of the gluon distribution.
     */
    const char* name();
    
    void write_pspace_grid(std::ostream& out);
    void write_rspace_grid(std::ostream& out);
    void write_satscale_grid(std::ostream& out);
protected:
    /**
     * Performs setup common to constructors
     */
    void setup(std::string pos_filename, std::string mom_filename, double xinit);
    /**
     * Initializes the variables used to compute Qs2(Y) from the position space function
     */
    void initialize_saturation_scale_from_position_space(double satscale_threshold);
    /**
     * Initializes the variables used to compute Qs2(Y) from the momentum space function
     */
    void initialize_saturation_scale_from_momentum_space(double satscale_threshold);
private:
    double r2min, r2max;
    double Yminr, Ymaxr;
    double q2min, q2max;
    double Yminp, Ymaxp;
    /** Values of ln(r2) for the interpolation. */
    double* r2_values;
    /** Values of Y for the position interpolation. */
    double* Y_values_rspace;
    /** Values of ln(q2) for the interpolation. */
    double* q2_values;
    /** Values of Y for the position interpolation. */
    double* Y_values_pspace;
    /** Values of the position space gluon distribution for interpolation. */
    double* S_dist;
    /** Values of the momentum space gluon distribution for interpolation. */
    double* F_dist;
    
    satscale_source_type satscale_source;
    
    // one of interp_dist_momentum_1D or interp_dist_momentum_2D will be NULL and the other one
    // will be used, depending on whether there is a range of Y values or just a single one
    gsl_interp* interp_dist_momentum_1D;
    interp2d* interp_dist_momentum_2D;
    // same for position
    gsl_interp* interp_dist_position_1D;
    interp2d* interp_dist_position_2D;
    
    // this will be used if we are extracting the saturation scale
    double* Qs2_values;
    gsl_interp* interp_Qs2_1D;
    
    gsl_interp_accel* r2_accel;
    gsl_interp_accel* q2_accel;
    gsl_interp_accel* Y_accel_r;
    gsl_interp_accel* Y_accel_p;
    
    size_t r2_dimension;
    size_t q2_dimension;
    size_t Y_dimension_r;
    size_t Y_dimension_p;
    
    size_t subinterval_limit;
    
    double Q02x0lambda;
    double lambda;
    
    std::string _name;
    
    friend class HybridGBWFileDataGluonDistribution;
    friend class HybridMVFileDataGluonDistribution;
};

class HybridGBWFileDataGluonDistribution : public FileDataGluonDistribution {
public:
    HybridGBWFileDataGluonDistribution(
        std::string pos_filename,
        std::string mom_filename,
        double Q02,
        double x0,
        double lambda,
        double xinit,
        enum satscale_source satscale_source,
        double satscale_threshold);
    HybridGBWFileDataGluonDistribution(
        std::string pos_filename,
        std::string mom_filename,
        double Q02,
        double x0,
        double lambda,
        double xinit);
    ~HybridGBWFileDataGluonDistribution();
    double S2(double r2, double Y);
    double S4(double r2, double s2, double t2, double Y);
    double Qs2(const double Y) const;
    double F(double q2, double Y);
private:
    GBWGluonDistribution gbw_dist;
};

class HybridMVFileDataGluonDistribution : public FileDataGluonDistribution {
public:
    HybridMVFileDataGluonDistribution(
        std::string pos_filename,
        std::string mom_filename,
        double LambdaMV,
        double gammaMV,
        double q2min,
        double q2max,
        double Ymin,
        double Ymax,
        double Q02,
        double x0,
        double lambda,
        double xinit,
        enum satscale_source satscale_source,
        double satscale_threshold,
        size_t subinterval_limit = 10000);
    HybridMVFileDataGluonDistribution(
        std::string pos_filename,
        std::string mom_filename,
        double LambdaMV,
        double gammaMV,
        double q2min,
        double q2max,
        double Ymin,
        double Ymax,
        double Q02,
        double x0,
        double lambda,
        double xinit,
        size_t subinterval_limit = 10000);
    ~HybridMVFileDataGluonDistribution();
    double S2(double r2, double Y);
    double S4(double r2, double s2, double t2, double Y);
    double Qs2(const double Y) const;
    double F(double q2, double Y);
private:
    MVGluonDistribution mv_dist;
};

/** Prints the name of the gluon distribution to the given output stream. */
std::ostream& operator<<(std::ostream& out, GluonDistribution& gdist);

/**
 * A wrapper object that prints out every call made to the gluon distribution.
 */
class GluonDistributionTraceWrapper : public GluonDistribution {
public:
    GluonDistributionTraceWrapper(GluonDistribution* gdist, const char* trace_filename = "trace_gdist.output") : GluonDistribution(), gdist(gdist), trace_stream(trace_filename) {}
    ~GluonDistributionTraceWrapper() { delete gdist; trace_stream.close(); }

    double S2(double r2, double Y);
    double S4(double r2, double s2, double t2, double Y);
    double F(double q2, double Y);
    double Qs2(const double Y) const;
    const char* name();
    
private:
    GluonDistribution* gdist;
    std::ofstream trace_stream;
};


#endif // _GLUONDIST_H_