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

#include <cassert>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_bessel.h>
#include "gluondist.h"
#include "interp2d.h"

using namespace std;

static NoGridException nogrid;

void GluonDistribution::write_pspace_grid(ostream& out) {
    throw nogrid;
}

void GluonDistribution::write_rspace_grid(ostream& out) {
    throw nogrid;
}

void GluonDistribution::write_satscale_grid(ostream& out) {
    throw nogrid;
}

GBWGluonDistribution::GBWGluonDistribution(double Q02, double x0, double lambda) : GluonDistribution(), Q02x0lambda(Q02 * pow(x0, lambda)), lambda(lambda) {
    ostringstream s;
    s << "GBW(Q02 = " << Q02 << ", x0 = " << x0 << ", lambda = " << lambda << ")";
    _name = s.str();
}

double GBWGluonDistribution::S2(double r2, double Y) {
    return exp(-0.25 * r2 * Qs2(Y));
}
double GBWGluonDistribution::S4(double r2, double s2, double t2, double Y) {
    return exp(-0.25 * Qs2(Y) * (s2 + t2));
}
double GBWGluonDistribution::F(double q2, double Y) {
    double _Qs2 = Qs2(Y);
    return M_1_PI * exp(-q2/_Qs2) / _Qs2;
}
double GBWGluonDistribution::Qs2(const double Y) const {
    return Q02x0lambda * exp(lambda * Y);
}
const char* GBWGluonDistribution::name() {
    return _name.c_str();
}


class PositionGDistIntegrationParameters {
public:
    double q; // use q instead of q^2 because that's what we need for the integrand
    double Qs2;
    size_t n; // for integrating a term of the series: this is the order of the term being integrated
    AbstractPositionGluonDistribution* gdist;
    PositionGDistIntegrationParameters(AbstractPositionGluonDistribution* gdist) : gdist(gdist), q(0), Qs2(0), n(0) {}
};

static double position_gdist_integrand(double r, void* closure) {
    PositionGDistIntegrationParameters* params = (PositionGDistIntegrationParameters*)closure;
    return 0.5 * M_1_PI * r * params->gdist->S2(r*r, params->Qs2) * gsl_sf_bessel_J0(params->q * r);
}

static double position_gdist_series_term_integrand(double r, void* closure) {
    PositionGDistIntegrationParameters* params = (PositionGDistIntegrationParameters*)closure;
    switch (params->n) {
        case 0:
            return 0.5 * M_1_PI * r * params->gdist->S2(r*r, params->Qs2);
        case 2:
            return -0.125 * M_1_PI * gsl_pow_3(r) * params->gdist->S2(r*r, params->Qs2);
        default: // a term not in the series
            return 0;
    }
}

AbstractPositionGluonDistribution::AbstractPositionGluonDistribution(double q2min, double q2max, double Ymin, double Ymax, size_t subinterval_limit) :
 GluonDistribution(),
 q2min(q2min), q2max(q2max), Ymin(Ymin), Ymax(Ymax),
 F_dist_leading_q2(NULL), F_dist_subleading_q2(NULL), F_dist(NULL),
 interp_dist_leading_q2(NULL), interp_dist_subleading_q2(NULL), interp_dist_momentum_1D(NULL), interp_dist_momentum_2D(NULL),
 q2_accel(NULL), Y_accel(NULL),
 q2_dimension(1), Y_dimension(1),
 subinterval_limit(subinterval_limit) {
 }
 
void AbstractPositionGluonDistribution::setup() {
    double step = 1.05;
    PositionGDistIntegrationParameters params(this);
    gsl_function func;
    func.params = &params;
    gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(subinterval_limit);
    
    double log_step = log(step);
    double log_q2min = log(q2min);
    double log_q2max = log(q2max);
    assert(log_step > 0);
    
    q2_dimension = (size_t)((log_q2max - log_q2min) / log_step) + 2; // subtracting logs rather than dividing may help accuracy
    while (q2_dimension < 4) { // 4 points needed for bicubic interpolation
        q2min /= step;
        log_q2min -= log_step;
        q2_dimension++;
    }
    if (Ymin != Ymax) {
        Y_dimension = (size_t)((Ymax - Ymin) / log_step) + 2;
        while (Y_dimension < 4) { // 4 points needed for bicubic interpolation
            Ymin -= log_step;
            Y_dimension++;
        }
    }
    
    log_q2_values = new double[q2_dimension];
    Y_values = new double[Y_dimension];
    
    double error; // throwaway
    
    // calculate the coefficients for the series approximation
    func.function = &position_gdist_series_term_integrand;
    F_dist_leading_q2 = new double[Y_dimension]; // zeroth order term in series around q2 = 0
    F_dist_subleading_q2 = new double[Y_dimension]; // second order term in series around q2 = 0
    for (size_t i_Y = 0; i_Y < Y_dimension; i_Y++) {
        Y_values[i_Y] = Ymin + i_Y * log_step;
        params.Qs2 = Qs2(Y_values[i_Y]); // TODO: verify whether sat scale is available at this point
        
        params.n = 0;
        gsl_integration_qagiu(&func, 0, 0, 0.0001, subinterval_limit, workspace, F_dist_leading_q2 + i_Y, &error);
        
        params.n = 2;
        gsl_integration_qagiu(&func, 0, 0, 0.0001, subinterval_limit, workspace, F_dist_subleading_q2 + i_Y, &error);
    }
    
    // calculate the values for the 2D interpolation
    func.function = &position_gdist_integrand;
    F_dist = new double[q2_dimension * Y_dimension];
    for (size_t i_q2 = 0; i_q2 < q2_dimension; i_q2++) {
        log_q2_values[i_q2] = log_q2min + i_q2 * log_step;
        params.q = exp(0.5 * log_q2_values[i_q2]);
        for (size_t i_Y = 0; i_Y < Y_dimension; i_Y++) {
            params.Qs2 = Qs2(Y_values[i_Y]);
            size_t index = INDEX_2D(i_q2, i_Y, q2_dimension, Qs2_dimension);
            gsl_integration_qagiu(&func, 0, 0, 0.0001, subinterval_limit, workspace, F_dist + index, &error);
        }
    }
    
    assert(log_q2_values[0] <= log_q2min);
    assert(log_q2_values[q2_dimension - 1] >= log_q2max);
    assert(Y_values[0] <= Ymin);
    assert(Y_values[Y_dimension - 1] >= Ymax);
    
    q2_accel = gsl_interp_accel_alloc();

    if (Y_dimension == 1) {
        interp_dist_momentum_1D = gsl_interp_alloc(gsl_interp_cspline, q2_dimension);
        gsl_interp_init(interp_dist_momentum_1D, log_q2_values, F_dist, q2_dimension);
    }
    else {
        Y_accel = gsl_interp_accel_alloc();
        
        interp_dist_leading_q2 = gsl_interp_alloc(gsl_interp_cspline, Y_dimension);
        gsl_interp_init(interp_dist_leading_q2, Y_values, F_dist_leading_q2, Y_dimension);
        
        interp_dist_subleading_q2 = gsl_interp_alloc(gsl_interp_cspline, Y_dimension);
        gsl_interp_init(interp_dist_subleading_q2, Y_values, F_dist_subleading_q2, Y_dimension);
        
        interp_dist_momentum_2D = interp2d_alloc(interp2d_bilinear, q2_dimension, Y_dimension);
        interp2d_init(interp_dist_momentum_2D, log_q2_values, Y_values, F_dist, q2_dimension, Y_dimension);
    }
}

AbstractPositionGluonDistribution::~AbstractPositionGluonDistribution() {
    delete[] F_dist;
    F_dist = NULL;
    delete[] F_dist_leading_q2;
    F_dist_leading_q2 = NULL;
    delete[] F_dist_subleading_q2;
    F_dist_subleading_q2 = NULL;
    interp2d_free(interp_dist_momentum_2D);
    interp_dist_momentum_2D = NULL;
    gsl_interp_free(interp_dist_leading_q2);
    interp_dist_momentum_2D = NULL;
    gsl_interp_free(interp_dist_subleading_q2);
    interp_dist_momentum_2D = NULL;
    gsl_interp_accel_free(q2_accel);
    q2_accel = NULL;
    gsl_interp_accel_free(Y_accel);
    Y_accel = NULL;
}

double AbstractPositionGluonDistribution::S4(double r2, double s2, double t2, double Y) {
    return S2(s2, Y) * S2(t2, Y);
}
double AbstractPositionGluonDistribution::F(double q2, double Y) {
    if (Y_dimension == 1) {
        if (q2 > q2min) {
            if (q2 > q2max || Y < Ymin || Y > Ymax) {
                throw GluonDistributionFRangeException(q2, Y);
            }
            return gsl_interp_eval(interp_dist_momentum_1D, log_q2_values, F_dist, log(q2), q2_accel);
        }
        else {
            double c0 = F_dist_leading_q2[0];
            double c2 = F_dist_subleading_q2[0];
            return c0 + c2 * q2;
        }
    }
    else {
        if (q2 > q2min) {
            return interp2d_eval(interp_dist_momentum_2D, log_q2_values, Y_values, F_dist, log(q2), Y, q2_accel, Y_accel);
        }
        else {
            double c0 = gsl_interp_eval(interp_dist_leading_q2, Y_values, F_dist_leading_q2, Y, Y_accel);
            double c2 = gsl_interp_eval(interp_dist_subleading_q2, Y_values, F_dist_subleading_q2, Y, Y_accel);
            return c0 + c2 * q2;
        }
    }
}

void AbstractPositionGluonDistribution::write_pspace_grid(ostream& out) {
    out << "q2\tQs2\tF" << endl;
    for (size_t i_q2 = 0; i_q2 < q2_dimension; i_q2++) {
        for (size_t i_Y = 0; i_Y < Y_dimension; i_Y++) {
            out << exp(log_q2_values[i_q2]) << "\t"
                << exp(Y_values[i_Y]) << "\t"
                << F_dist[INDEX_2D(i_q2, i_Y, q2_dimension, Y_dimension)] << endl;
        }
    }
}

MVGluonDistribution::MVGluonDistribution(double LambdaMV, double gammaMV, double q2min, double q2max, double Ymin, double Ymax, double Q02, double x0, double lambda, size_t subinterval_limit) :
  AbstractPositionGluonDistribution(q2min, q2max, Ymin, Ymax, subinterval_limit), LambdaMV(LambdaMV), gammaMV(gammaMV), Q02x0lambda(Q02 * pow(x0, lambda)), lambda(lambda) {
    ostringstream s;
    s << "MV(LambdaMV = " << LambdaMV << ", gammaMV = " << gammaMV << ", q2min = " << q2min << ", q2max = " << q2max << ", Ymin = " << Ymin << ", Ymax = " << Ymax << ", Q02 = " << Q02 << ", x0 = " << x0 << ", lambda = " << lambda << ")";
    _name = s.str();
    setup();
}

double MVGluonDistribution::S2(double r2, double Y) {
    return pow(M_E + 1.0 / (sqrt(r2) * LambdaMV), -0.25 * pow(r2 * Qs2(Y), gammaMV));
}

double MVGluonDistribution::Qs2(const double Y) const {
    return Q02x0lambda * exp(lambda * Y);
}

const char* MVGluonDistribution::name() {
    return _name.c_str();
}

FixedSaturationMVGluonDistribution::FixedSaturationMVGluonDistribution(double LambdaMV, double gammaMV, double q2min, double q2max, double YMV, double Q02, double x0, double lambda, size_t subinterval_limit) :
  MVGluonDistribution(LambdaMV, gammaMV, q2min, q2max, YMV, YMV, Q02, x0, lambda, subinterval_limit) {
    // calculate_interpolation_grid runs in the superclass constructor
    ostringstream s;
    s << "fMV(LambdaMV = " << LambdaMV << ", gammaMV = " << gammaMV << ", q2min = " << q2min << ", q2max = " << q2max << ", YMV = " << YMV << ", Q02 = " << Q02 << ", x0 = " << x0 << ", lambda = " << lambda << ")";
    _name = s.str();
}

double FixedSaturationMVGluonDistribution::S2(double r2, double Y) {
    return MVGluonDistribution::S2(r2, YMV);
}

struct double_triplet {
    double x, y, z;
    bool processed;
};

bool index_comparator(const double_triplet& a, const double_triplet& b) {
    // be generic in case it changes
    const static size_t row_major = 1, col_major = 2;
    // row index <-> x index
    // col index <-> y index
    switch (INDEX_2D(row_major, col_major, 0, 0)) {
        case row_major:
            return a.y < b.y || (a.y == b.y && a.x < b.x);
        case col_major:
            return a.x < b.x || (a.x == b.x && a.y < b.y);
    }
}

void read_from_file(const string filename, size_t& x_dimension, size_t& y_dimension, double*& x_values, double*& y_values, double*& z_values) {
    vector<double_triplet> file_data;
    set<double> xvals, yvals;

    ifstream in(filename.c_str());
    double_triplet point;
    point.processed = false;
    while (in >> point.x >> point.y >> point.z) {
        if (gsl_isnan(point.z)) {
            GSL_ERROR_VOID("NaN in data file", GSL_EINVAL);
        }
        file_data.push_back(point);
        xvals.insert(point.x);
        yvals.insert(point.y);
    }
    in.close();
    
    x_dimension = xvals.size();
    y_dimension = yvals.size();
    if (x_dimension * y_dimension != file_data.size()) {
        GSL_ERROR_VOID("Points not on a grid", GSL_EINVAL);
    }
    x_values = new double[x_dimension];
    y_values = new double[y_dimension];
    z_values = new double[x_dimension * y_dimension];
    
    sort(file_data.begin(), file_data.end(), index_comparator);
    
    copy(xvals.begin(), xvals.end(), x_values);
    copy(yvals.begin(), yvals.end(), y_values);
    for (size_t ix = 0; ix < x_dimension; ix++) {
        for (size_t iy = 0; iy < y_dimension; iy++) {
            size_t index = INDEX_2D(ix, iy, x_dimension, y_dimension);
            if (file_data[index].processed) {
                GSL_ERROR_VOID("Duplicate x, y in file", GSL_EINVAL);
            }
            if (x_values[ix] != file_data[index].x) {
                GSL_ERROR_VOID("Points out of order in gdist file", GSL_EINVAL);
            }
            if (y_values[iy] != file_data[index].y) {
                GSL_ERROR_VOID("Points out of order in gdist file", GSL_EINVAL);
            }
            z_values[index] = file_data[index].z;
            file_data[index].processed = true;
        }
    }
}

FileDataGluonDistribution::FileDataGluonDistribution(string pos_filename, string mom_filename, double xinit, enum satscale_source satscale_source, double satscale_threshold) :
  GluonDistribution(), Q02x0lambda(0), lambda(0), Qs2_values(NULL), satscale_source(satscale_source) {
    setup(pos_filename, mom_filename, xinit);
    ostringstream s;
    s << "file(pos_filename = " << pos_filename << ", mom_filename = " << mom_filename << ", xinit = " << xinit;
    switch (satscale_source) {
        case MOMENTUM_THRESHOLD:
            initialize_saturation_scale_from_momentum_space(satscale_threshold);
            s << ", extracted saturation scale from momentum)";
            break;
        case POSITION_THRESHOLD:
            initialize_saturation_scale_from_position_space(satscale_threshold);
            s << ", extracted saturation scale from position)";
            break;
        default:
            throw "Invalid saturation scale source";
    }
    _name = s.str();
}


FileDataGluonDistribution::FileDataGluonDistribution(string pos_filename, string mom_filename, double Q02, double x0, double lambda, double xinit) :
  GluonDistribution(), Q02x0lambda(Q02 * pow(x0, lambda)), lambda(lambda), Qs2_values(NULL), satscale_source(NONE) {
    setup(pos_filename, mom_filename, xinit);
    ostringstream s;
    s << "file(pos_filename = " << pos_filename << ", mom_filename = " << mom_filename << ", xinit = " << xinit << ", default saturation scale)";
    _name = s.str();
}

FileDataGluonDistribution::~FileDataGluonDistribution() {
    delete[] r2_values, q2_values, Y_values_rspace, Y_values_pspace, S_dist, F_dist, Qs2_values;
    gsl_interp_accel_free(r2_accel);
    gsl_interp_accel_free(q2_accel);
    gsl_interp_accel_free(Y_accel_r);
    gsl_interp_accel_free(Y_accel_p);
    if (Y_dimension_r == 1) {
        gsl_interp_free(interp_dist_position_1D);
        gsl_interp_free(interp_dist_momentum_1D);
    }
    else {
        interp2d_free(interp_dist_position_2D);
        interp2d_free(interp_dist_momentum_2D);
        gsl_interp_free(interp_Qs2_1D);
    }
}

void FileDataGluonDistribution::setup(string pos_filename, string mom_filename, double xinit) {
    // the rapidity values we read from the file are considered values of ΔY
    // relative to some initial rapidity, Yinit.
    read_from_file(pos_filename, r2_dimension, Y_dimension_r, r2_values, Y_values_rspace, S_dist);
    read_from_file(mom_filename, q2_dimension, Y_dimension_p, q2_values, Y_values_pspace, F_dist);
    
    // apply the offset
    if (xinit != 1) {
        double Yinit = -log(xinit);
        for (size_t i = 0; i < Y_dimension_r; i++) {
            Y_values_rspace[i] += Yinit;
        }
        for (size_t i = 0; i < Y_dimension_p; i++) {
            Y_values_pspace[i] += Yinit;
        }
    }
    
    r2min = r2_values[0];
    r2max = r2_values[r2_dimension-1];
    Yminr = Y_values_rspace[0];
    Ymaxr = Y_values_rspace[Y_dimension_r-1];
    q2min = q2_values[0];
    q2max = q2_values[q2_dimension-1];
    Yminp = Y_values_pspace[0];
    Ymaxp = Y_values_pspace[Y_dimension_p-1];
    
    r2_accel = gsl_interp_accel_alloc();
    q2_accel = gsl_interp_accel_alloc();

    if (Y_dimension_r == 1) {
        assert(Y_dimension_p == 1);
        assert(Yminr == Ymaxr);
        assert(Yminp == Ymaxp);
        interp_dist_position_1D = gsl_interp_alloc(gsl_interp_cspline, r2_dimension);
        gsl_interp_init(interp_dist_position_1D, r2_values, S_dist, r2_dimension);
        interp_dist_momentum_1D = gsl_interp_alloc(gsl_interp_cspline, q2_dimension);
        gsl_interp_init(interp_dist_momentum_1D, q2_values, F_dist, q2_dimension);
    }
    else {
        assert(Y_dimension_p > 1);
        Y_accel_r = gsl_interp_accel_alloc();
        Y_accel_p = gsl_interp_accel_alloc();

        interp_dist_position_2D = interp2d_alloc(interp2d_bilinear, r2_dimension, Y_dimension_r);
        interp2d_init(interp_dist_position_2D, r2_values, Y_values_rspace, S_dist, r2_dimension, Y_dimension_r);
        
        interp_dist_momentum_2D = interp2d_alloc(interp2d_bilinear, q2_dimension, Y_dimension_p);
        interp2d_init(interp_dist_momentum_2D, q2_values, Y_values_pspace, F_dist, q2_dimension, Y_dimension_p);
    }
}

struct EvaluationParameters {
    FileDataGluonDistribution* gdist;
    double Y;
    double threshold;
};

double bracket_root(double (*f)(double, void*), EvaluationParameters* params, double min, double max) {
    const static size_t max_iter = 1000; // should be way more than necessary
    size_t iter = 0;
    int status;
    double current_root;
    
    { // check that the function has opposite signs at the endpoints
        double minv = f(min, params);
        double maxv = f(max, params);
        if (minv == 0) {
            return min;
        }
        if (maxv == 0) {
            return max;
        }
        if ((minv > 0 && maxv > 0) || (minv < 0 && maxv < 0)) {
            GSL_ERROR_VAL("Zero crossing not found", GSL_EINVAL, 0.0);
        }
    }
    
    gsl_function func;
    func.function = f;
    func.params = params;
    
    gsl_root_fsolver* solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    gsl_root_fsolver_set(solver, &func, min, max);
    
    do {
        status = gsl_root_fsolver_iterate(solver);
        current_root = gsl_root_fsolver_root(solver);
        min = gsl_root_fsolver_x_lower(solver);
        max = gsl_root_fsolver_x_upper(solver);
        status = gsl_root_test_interval(min, max, 0, 1e-3);
    } while (status == GSL_CONTINUE && ++iter < max_iter);
    gsl_root_fsolver_free(solver);
    
    return current_root;
}

double evaluate_rspace_threshold_criterion(double r2, void* params) {
    EvaluationParameters* p = (EvaluationParameters*)params;
    return p->gdist->S2(r2, p->Y) - p->threshold;
}

void FileDataGluonDistribution::initialize_saturation_scale_from_position_space(double satscale_threshold) {
    // strategy: search for the point where S(r,Y) = T
    assert(Y_dimension_r >= 1);
    
    double last_Rs2;
    EvaluationParameters p;
    p.gdist = this;
    p.threshold = satscale_threshold;
    
    Qs2_values = new double[Y_dimension_r];
    for (size_t i = 0; i < Y_dimension_r; i++) {
        p.Y = Y_values_rspace[i];
        last_Rs2 = bracket_root(evaluate_rspace_threshold_criterion, &p, r2min, r2max);
        Qs2_values[i] = 1.0 / last_Rs2;
    }
    if (Y_dimension_r > 1) {
        interp_Qs2_1D = gsl_interp_alloc(gsl_interp_cspline, Y_dimension_r);
        gsl_interp_init(interp_Qs2_1D, Y_values_rspace, Qs2_values, Y_dimension_r);
    }
}

double evaluate_pspace_threshold_criterion(double k2, void* params) {
    EvaluationParameters* p = (EvaluationParameters*)params;
    return p->gdist->F(k2, p->Y) - p->threshold;
}

void FileDataGluonDistribution::initialize_saturation_scale_from_momentum_space(double satscale_threshold) {
    // strategy: search for the point where F(k,Y) = T F(k0, Y) 
    assert(Y_dimension_p >= 1);
    
    double last_Qs2;
    EvaluationParameters p;
    p.gdist = this;
    
    Qs2_values = new double[Y_dimension_p];
    for (size_t i = 0; i < Y_dimension_p; i++) {
        p.Y = Y_values_pspace[i];
        p.threshold = satscale_threshold * F(q2min, p.Y);
        Qs2_values[i] = last_Qs2 = bracket_root(evaluate_pspace_threshold_criterion, &p, q2min, q2max);
    }
    if (Y_dimension_p > 1) {
        interp_Qs2_1D = gsl_interp_alloc(gsl_interp_cspline, Y_dimension_p);
        gsl_interp_init(interp_Qs2_1D, Y_values_pspace, Qs2_values, Y_dimension_p);
    }
}

double FileDataGluonDistribution::S2(double r2, double Y) {
    if (r2 < r2min || r2 > r2max || Y < Yminr || Y > Ymaxr) {
        throw GluonDistributionS2RangeException(r2, Y);
    }
    if (Y_dimension_r == 1) {
        return gsl_interp_eval(interp_dist_position_1D, r2_values, S_dist, r2, r2_accel);
    }
    else {
        return interp2d_eval(interp_dist_position_2D, r2_values, Y_values_rspace, S_dist, r2, Y, r2_accel, Y_accel_r);
    }
}

double FileDataGluonDistribution::S4(double r2, double s2, double t2, double Y) {
    if (s2 < r2min || s2 > r2max || t2 < r2min || t2 > r2max || Y < Yminr || Y > Ymaxr) {
        throw GluonDistributionS4RangeException(r2, s2, t2, Y);
    }
    return S2(s2, Y) * S2(t2, Y);
}

double FileDataGluonDistribution::F(double q2, double Y) {
    if (q2 < q2min || q2 > q2max || Y < Yminr || Y > Ymaxr) {
        throw GluonDistributionFRangeException(q2, Y);
    }
    if (Y_dimension_p == 1) {
        return gsl_interp_eval(interp_dist_momentum_1D, q2_values, F_dist, q2, q2_accel);
    }
    else {
        return interp2d_eval(interp_dist_momentum_2D, q2_values, Y_values_pspace, F_dist, q2, Y, q2_accel, Y_accel_p);
    }
}

double FileDataGluonDistribution::Fprime(double q2, double Y) {
    if (q2 < q2min || q2 > q2max || Y < Yminr || Y > Ymaxr) {
        throw GluonDistributionFRangeException(q2, Y);
    }
    if (Y_dimension_p == 1) {
        return gsl_interp_eval_deriv(interp_dist_momentum_1D, q2_values, F_dist, q2, q2_accel);
    }
    else {
        return interp2d_eval_deriv_x(interp_dist_momentum_2D, q2_values, Y_values_pspace, F_dist, q2, Y, q2_accel, Y_accel_p);
    }
}

double FileDataGluonDistribution::Fpprime(double q2, double Y) {
    if (q2 < q2min || q2 > q2max || Y < Yminr || Y > Ymaxr) {
        throw GluonDistributionFRangeException(q2, Y);
    }
    if (Y_dimension_p == 1) {
        return gsl_interp_eval_deriv2(interp_dist_momentum_1D, q2_values, F_dist, q2, q2_accel);
    }
    else {
        return interp2d_eval_deriv_xx(interp_dist_momentum_2D, q2_values, Y_values_pspace, F_dist, q2, Y, q2_accel, Y_accel_p);
    }
}

double FileDataGluonDistribution::S2prime(double r2, double Y) {
    if (r2 < r2min || r2 > r2max || Y < Yminr || Y > Ymaxr) {
        throw GluonDistributionS2RangeException(r2, Y);
    }
    if (Y_dimension_r == 1) {
        return gsl_interp_eval_deriv(interp_dist_position_1D, r2_values, S_dist, r2, r2_accel);
    }
    else {
        return interp2d_eval_deriv_x(interp_dist_position_2D, r2_values, Y_values_rspace, S_dist, r2, Y, r2_accel, Y_accel_r);
    }
}

double FileDataGluonDistribution::S2pprime(double r2, double Y) {
    if (r2 < r2min || r2 > r2max || Y < Yminr || Y > Ymaxr) {
        throw GluonDistributionS2RangeException(r2, Y);
    }
    if (Y_dimension_r == 1) {
        return gsl_interp_eval_deriv2(interp_dist_position_1D, r2_values, S_dist, r2, r2_accel);
    }
    else {
        return interp2d_eval_deriv_xx(interp_dist_position_2D, r2_values, Y_values_rspace, S_dist, r2, Y, r2_accel, Y_accel_r);
    }
}

double FileDataGluonDistribution::Qs2(const double Y) const {
    if (Qs2_values == NULL) {
        assert(Q02x0lambda > 0);
        assert(lambda > 0);
        return Q02x0lambda * exp(lambda * Y);
    }
    else if (Y_dimension_p == 1) {
        assert(Y_dimension_r == 1);
        return Qs2_values[0];
    }
    else {
        switch (satscale_source) {
            case MOMENTUM_THRESHOLD:
                return gsl_interp_eval(interp_Qs2_1D, Y_values_pspace, Qs2_values, Y, Y_accel_p);
            case POSITION_THRESHOLD:
                return gsl_interp_eval(interp_Qs2_1D, Y_values_rspace, Qs2_values, Y, Y_accel_r);
            default:
                assert(false);
        }
    }
}

const char* FileDataGluonDistribution::name() {
    return _name.c_str();
}

void FileDataGluonDistribution::write_pspace_grid(ostream& out) {
    out << "q2\tQs2\tF" << endl;
    for (size_t i_q2 = 0; i_q2 < q2_dimension; i_q2++) {
        for (size_t i_Y = 0; i_Y < Y_dimension_p; i_Y++) {
            out << q2_values[i_q2] << "\t"
                << Y_values_pspace[i_Y] << "\t"
                << F_dist[INDEX_2D(i_q2, i_Y, q2_dimension, Y_dimension_p)] << endl;
        }
    }
}

void FileDataGluonDistribution::write_rspace_grid(ostream& out) {
    out << "r2\tQs2\tS" << endl;
    for (size_t i_r2 = 0; i_r2 < r2_dimension; i_r2++) {
        for (size_t i_Y = 0; i_Y < Y_dimension_r; i_Y++) {
            out << r2_values[i_r2] << "\t"
                << Y_values_rspace[i_Y] << "\t"
                << F_dist[INDEX_2D(i_r2, i_Y, r2_dimension, Y_dimension_r)] << endl;
        }
    }
}

void FileDataGluonDistribution::write_satscale_grid(ostream& out) {
    if (Qs2_values == NULL) {
        throw nogrid;
    }
    size_t Y_dimension;
    double* Y_values;
    switch (satscale_source) {
        case MOMENTUM_THRESHOLD:
            Y_dimension = Y_dimension_p;
            Y_values = Y_values_pspace;
            break;
        case POSITION_THRESHOLD:
            Y_dimension = Y_dimension_r;
            Y_values = Y_values_rspace;
            break;
    }
    out << "Y\tx\tQs2" << endl;
    for (size_t i_Y = 0; i_Y < Y_dimension; i_Y++) {
        out << Y_values[i_Y] << "\t"
            << exp(-Y_values[i_Y]) << "\t"
            << Qs2_values[i_Y] << endl;
    }
}

HybridGBWFileDataGluonDistribution::HybridGBWFileDataGluonDistribution(
 string pos_filename,
 string mom_filename,
 double Q02, double x0, double lambda,
 double xinit,
 enum satscale_source satscale_source,
 double satscale_threshold) :
  FileDataGluonDistribution(pos_filename, mom_filename, xinit, satscale_source, satscale_threshold),
  gbw_dist(Q02, x0, lambda) {
}

HybridGBWFileDataGluonDistribution::HybridGBWFileDataGluonDistribution(
 string pos_filename,
 string mom_filename,
 double Q02, double x0, double lambda,
 double xinit) :
  FileDataGluonDistribution(pos_filename, mom_filename, Q02, x0, lambda, xinit),
  gbw_dist(Q02, x0, lambda) {
}

HybridGBWFileDataGluonDistribution::~HybridGBWFileDataGluonDistribution() {}

double HybridGBWFileDataGluonDistribution::S2(double r2, double Y) {
    if (r2 < r2min) {
        return gbw_dist.S2(r2, Y);
    }
    else {
        return FileDataGluonDistribution::S2(r2, Y);
    }
}

double HybridGBWFileDataGluonDistribution::F(double q2, double Y) {
    if (q2 < q2min) {
        return gbw_dist.F(q2, Y);
    }
    else {
        return FileDataGluonDistribution::F(q2, Y);
    }
}

double HybridGBWFileDataGluonDistribution::Qs2(const double Y) const {
    if (satscale_source == POSITION_THRESHOLD && Y < Yminr) {
        return gbw_dist.Qs2(Y);
    }
    else if (satscale_source == MOMENTUM_THRESHOLD && Y < Yminp) {
        return gbw_dist.Qs2(Y);
    }
    return FileDataGluonDistribution::Qs2(Y);
}

ostream& operator<<(ostream& out, GluonDistribution& gdist) {
    out << gdist.name();
    return out;
}

double GluonDistributionTraceWrapper::F(double q2, double Y) {
    double val = gdist->F(q2, Y);
    trace_stream << "F\t" << q2 << "\t" << Y << "\t" << val << endl;
    return val;
}

double GluonDistributionTraceWrapper::S2(double r2, double Y) {
    double val = gdist->S2(r2, Y);
    trace_stream << "S2\t" << r2 << "\t" << Y << "\t" << val << endl;
    return val;
}

double GluonDistributionTraceWrapper::S4(double r2, double s2, double t2, double Y) {
    double val = gdist->S4(r2, s2, t2, Y);
    trace_stream << "S4\t" << r2 << "\t" << s2 << "\t" << t2 << "\t" << Y << "\t" << val << endl;
    return val;
}

double GluonDistributionTraceWrapper::Qs2(const double Y) const {
    return gdist->Qs2(Y);
}

const char* GluonDistributionTraceWrapper::name() {
    return gdist->name();
}
