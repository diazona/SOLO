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
#include <gsl/gsl_sf_bessel.h>
#include "gluondist.h"
#include "interp2d.h"

using namespace std;

double GBWGluonDistribution::S2(double r2, double Qs2) {
    return exp(-0.25 * r2 * Qs2);
}
double GBWGluonDistribution::S4(double r2, double s2, double t2, double Qs2) {
    return exp(-0.25 * Qs2 * (s2 + t2));
}
double GBWGluonDistribution::F(double q2, double Qs2) {
    return M_1_PI * exp(-q2/Qs2) / Qs2;
}
const char* GBWGluonDistribution::name() {
    return "GBW";
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

AbstractPositionGluonDistribution::AbstractPositionGluonDistribution(double q2min, double q2max, double Qs2min, double Qs2max, size_t subinterval_limit) :
 q2min(q2min), q2max(q2max), Qs2min(Qs2min), Qs2max(Qs2max),
 F_dist_leading_q2(NULL), F_dist_subleading_q2(NULL), F_dist(NULL),
 interp_dist_leading_q2(NULL), interp_dist_subleading_q2(NULL), interp_dist_momentum_1D(NULL), interp_dist_momentum_2D(NULL),
 q2_accel(NULL), Qs2_accel(NULL),
 q2_dimension(1), Qs2_dimension(1),
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
    double log_Qs2min = log(Qs2min);
    double log_q2max = log(q2max);
    double log_Qs2max = log(Qs2max);
    assert(log_step > 0);
    
    q2_dimension = (size_t)((log_q2max - log_q2min) / log_step) + 2; // subtracting logs rather than dividing may help accuracy
    while (q2_dimension < 4) { // 4 points needed for bicubic interpolation
        q2min /= step;
        log_q2min -= log_step;
        q2_dimension++;
    }
    if (Qs2min != Qs2max) {
        Qs2_dimension = (size_t)((log_Qs2max - log_Qs2min) / log_step) + 2;
        while (Qs2_dimension < 4) { // 4 points needed for bicubic interpolation
            Qs2min /= step;
            log_Qs2min -= log_step;
            Qs2_dimension++;
        }
    }
    
    log_q2_values = new double[q2_dimension];
    log_Qs2_values = new double[Qs2_dimension];
    
    double error; // throwaway
    
    // calculate the coefficients for the series approximation
    func.function = &position_gdist_series_term_integrand;
    F_dist_leading_q2 = new double[Qs2_dimension]; // zeroth order term in series around q2 = 0
    F_dist_subleading_q2 = new double[Qs2_dimension]; // second order term in series around q2 = 0
    for (size_t i_Qs2 = 0; i_Qs2 < Qs2_dimension; i_Qs2++) {
        log_Qs2_values[i_Qs2] = log_Qs2min + i_Qs2 * log_step;
        params.Qs2 = exp(log_Qs2_values[i_Qs2]);
        
        params.n = 0;
        gsl_integration_qagiu(&func, 0, 0, 0.0001, subinterval_limit, workspace, F_dist_leading_q2 + i_Qs2, &error);
        
        params.n = 2;
        gsl_integration_qagiu(&func, 0, 0, 0.0001, subinterval_limit, workspace, F_dist_subleading_q2 + i_Qs2, &error);
    }
    
    // calculate the values for the 2D interpolation
    func.function = &position_gdist_integrand;
    F_dist = new double[q2_dimension * Qs2_dimension];
    for (size_t i_q2 = 0; i_q2 < q2_dimension; i_q2++) {
        log_q2_values[i_q2] = log_q2min + i_q2 * log_step;
        params.q = exp(0.5 * log_q2_values[i_q2]);
        for (size_t i_Qs2 = 0; i_Qs2 < Qs2_dimension; i_Qs2++) {
            params.Qs2 = exp(log_Qs2_values[i_Qs2]);
            size_t index = INDEX_2D(i_q2, i_Qs2, q2_dimension, Qs2_dimension);
            gsl_integration_qagiu(&func, 0, 0, 0.0001, subinterval_limit, workspace, F_dist + index, &error);
        }
    }
    
    assert(log_q2_values[0] <= log_q2min);
    assert(log_q2_values[q2_dimension - 1] >= log_q2max);
    assert(log_Qs2_values[0] <= log_Qs2min);
    assert(log_Qs2_values[Qs2_dimension - 1] >= log_Qs2max);
    
    q2_accel = gsl_interp_accel_alloc();

    if (Qs2_dimension == 1) {
        interp_dist_momentum_1D = gsl_interp_alloc(gsl_interp_cspline, q2_dimension);
        gsl_interp_init(interp_dist_momentum_1D, log_q2_values, F_dist, q2_dimension);
    }
    else {
        Qs2_accel = gsl_interp_accel_alloc();
        
        interp_dist_leading_q2 = gsl_interp_alloc(gsl_interp_cspline, Qs2_dimension);
        gsl_interp_init(interp_dist_leading_q2, log_Qs2_values, F_dist_leading_q2, Qs2_dimension);
        
        interp_dist_subleading_q2 = gsl_interp_alloc(gsl_interp_cspline, Qs2_dimension);
        gsl_interp_init(interp_dist_subleading_q2, log_Qs2_values, F_dist_subleading_q2, Qs2_dimension);
        
        interp_dist_momentum_2D = interp2d_alloc(interp2d_bilinear, q2_dimension, Qs2_dimension);
        interp2d_init(interp_dist_momentum_2D, log_q2_values, log_Qs2_values, F_dist, q2_dimension, Qs2_dimension);
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
    gsl_interp_accel_free(Qs2_accel);
    Qs2_accel = NULL;
}

double AbstractPositionGluonDistribution::S4(double r2, double s2, double t2, double Qs2) {
    return S2(s2, Qs2) * S2(t2, Qs2);
}
double AbstractPositionGluonDistribution::F(double q2, double Qs2) {
    if (Qs2_dimension == 1) {
        if (q2 > q2min) {
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
            return interp2d_eval(interp_dist_momentum_2D, log_q2_values, log_Qs2_values, F_dist, log(q2), log(Qs2), q2_accel, Qs2_accel);
        }
        else {
            double c0 = gsl_interp_eval(interp_dist_leading_q2, log_Qs2_values, F_dist_leading_q2, log(Qs2), Qs2_accel);
            double c2 = gsl_interp_eval(interp_dist_subleading_q2, log_Qs2_values, F_dist_subleading_q2, log(Qs2), Qs2_accel);
            return c0 + c2 * q2;
        }
    }
}

MVGluonDistribution::MVGluonDistribution(double LambdaMV, double gammaMV, double q2min, double q2max, double Qs2min, double Qs2max, size_t subinterval_limit) :
  AbstractPositionGluonDistribution(q2min, q2max, Qs2min, Qs2max, subinterval_limit), LambdaMV(LambdaMV), gammaMV(gammaMV) {
    ostringstream s;
    s << "MV(LambdaMV = " << LambdaMV << ", gammaMV = " << gammaMV << ", q2min = " << q2min << ", q2max = " << q2max << ", Qs2min = " << Qs2min << ", Qs2max = " << Qs2max << ")";
    _name = s.str();
    setup();
}

double MVGluonDistribution::S2(double r2, double Qs2) {
    return pow(M_E + 1.0 / (sqrt(r2) * LambdaMV), -0.25 * pow(r2 * Qs2, gammaMV));
}

const char* MVGluonDistribution::name() {
    return _name.c_str();
}

FixedSaturationMVGluonDistribution::FixedSaturationMVGluonDistribution(double LambdaMV, double gammaMV, double q2min, double q2max, double Qs02MV, size_t subinterval_limit) :
  MVGluonDistribution(LambdaMV, gammaMV, q2min, q2max, Qs02MV, Qs02MV, subinterval_limit) {
    // calculate_interpolation_grid runs in the superclass constructor
    ostringstream s;
    s << "fMV(LambdaMV = " << LambdaMV << ", gammaMV = " << gammaMV << ", q2min = " << q2min << ", q2max = " << q2max << ", Qs02MV = " << Qs02MV << ")";
    _name = s.str();
}

double FixedSaturationMVGluonDistribution::S2(double r2, double Qs2) {
    return MVGluonDistribution::S2(r2, Qs02MV);
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

FileDataGluonDistribution::FileDataGluonDistribution(string pos_filename, string mom_filename, double Q02, double x0, double lambda) {
    double Q02x0lambda = Q02 * pow(x0, lambda);
    double* Y_values_rspace;
    double* Y_values_pspace;
    read_from_file(pos_filename, r2_dimension, Qs2_dimension_r, r2_values, Y_values_rspace, S_dist);
    read_from_file(mom_filename, q2_dimension, Qs2_dimension_p, q2_values, Y_values_pspace, F_dist);
    Qs2_values_rspace = new double[Qs2_dimension_r];
    Qs2_values_pspace = new double[Qs2_dimension_p];
    for (size_t i = 0; i < Qs2_dimension_r; i++) {
        Qs2_values_rspace[i] = Q02x0lambda * exp(lambda * Y_values_rspace[i]);
    }
    for (size_t i = 0; i < Qs2_dimension_p; i++) {
        Qs2_values_pspace[i] = Q02x0lambda * exp(lambda * Y_values_pspace[i]);
    }
    delete[] Y_values_rspace;
    Y_values_rspace = NULL;
    delete[] Y_values_pspace;
    Y_values_pspace = NULL;
    
    r2_accel = gsl_interp_accel_alloc();
    q2_accel = gsl_interp_accel_alloc();

    if (Qs2_dimension_r == 1) {
        assert(Qs2_dimension_p == 1);
        interp_dist_position_1D = gsl_interp_alloc(gsl_interp_cspline, r2_dimension);
        gsl_interp_init(interp_dist_position_1D, r2_values, S_dist, r2_dimension);
        interp_dist_momentum_1D = gsl_interp_alloc(gsl_interp_cspline, q2_dimension);
        gsl_interp_init(interp_dist_momentum_1D, q2_values, S_dist, q2_dimension);
    }
    else {
        assert(Qs2_dimension_p > 1);
        Qs2_accel_r = gsl_interp_accel_alloc();
        Qs2_accel_p = gsl_interp_accel_alloc();

        interp_dist_position_2D = interp2d_alloc(interp2d_bilinear, r2_dimension, Qs2_dimension_r);
        interp2d_init(interp_dist_position_2D, r2_values, Qs2_values_rspace, S_dist, r2_dimension, Qs2_dimension_r);
        
        interp_dist_momentum_2D = interp2d_alloc(interp2d_bilinear, q2_dimension, Qs2_dimension_p);
        interp2d_init(interp_dist_momentum_2D, q2_values, Qs2_values_pspace, F_dist, q2_dimension, Qs2_dimension_p);
    }
    
    ostringstream s;
    s << "file(pos_filename = " << pos_filename << ", mom_filename = " << mom_filename << ")";
    _name = s.str();
}

FileDataGluonDistribution::~FileDataGluonDistribution() {
    delete[] r2_values, q2_values, Qs2_values_rspace, Qs2_values_pspace, S_dist, F_dist;
    gsl_interp_accel_free(r2_accel);
    gsl_interp_accel_free(q2_accel);
    gsl_interp_accel_free(Qs2_accel_r);
    gsl_interp_accel_free(Qs2_accel_p);
    if (Qs2_dimension_r == 1) {
        gsl_interp_free(interp_dist_position_1D);
        gsl_interp_free(interp_dist_momentum_1D);
    }
    else {
        interp2d_free(interp_dist_position_2D);
        interp2d_free(interp_dist_momentum_2D);
    }
}

double FileDataGluonDistribution::S2(double r2, double Qs2) {
    if (Qs2_dimension_r == 1) {
        return gsl_interp_eval(interp_dist_position_1D, r2_values, F_dist, r2, r2_accel);
    }
    else {
        return interp2d_eval(interp_dist_position_2D, r2_values, Qs2_values_rspace, F_dist, r2, Qs2, r2_accel, Qs2_accel_r);
    }
}

double FileDataGluonDistribution::S4(double r2, double s2, double t2, double Qs2) {
    return S2(s2, Qs2) * S2(t2, Qs2);
}

double FileDataGluonDistribution::F(double q2, double Qs2) {
    if (Qs2_dimension_p == 1) {
        return gsl_interp_eval(interp_dist_momentum_1D, q2_values, F_dist, q2, q2_accel);
    }
    else {
        return interp2d_eval(interp_dist_momentum_2D, q2_values, Qs2_values_pspace, F_dist, q2, Qs2, q2_accel, Qs2_accel_p);
    }
}

const char* FileDataGluonDistribution::name() {
    return _name.c_str();
}


ostream& operator<<(ostream& out, GluonDistribution& gdist) {
    out << gdist.name();
    return out;
}

#ifdef GLUON_DIST_DRIVER
#include <cstdlib>
#include <iostream>

void AbstractPositionGluonDistribution::write_grid() {
    cout << "q2\tQs2\tF" << endl;
    for (size_t i_q2 = 0; i_q2 < q2_dimension; i_q2++) {
        for (size_t i_Qs2 = 0; i_Qs2 < Qs2_dimension; i_Qs2++) {
            cout << exp(log_q2_values[i_q2]) << "\t"
                    << exp(log_Qs2_values[i_Qs2]) << "\t"
                    << F_dist[INDEX_2D(i_q2, i_Qs2, q2_dimension, Qs2_dimension)] << endl;
        }
    }
}

void handle_input(AbstractPositionGluonDistribution& gdist) {
    cin.peek();
    if (cin.eof()) {
        // write out the grid
        gdist.write_grid();
    }
    else {
        double q2, Qs2;
        cin >> q2 >> Qs2;
        while (!cin.eof()) {
            cout << q2 << "\t"<< Qs2 << "\t" << gdist.F(q2, Qs2) << endl;
            cin >> q2 >> Qs2;
        }
    }
}

/**
 * A driver program that constructs a gluon distribution object
 * and prints its grid to standard output.
 */
int main(int argc, char** argv) {
    if (argc < 8) {
        cerr << "Needs 7 arguments: gdist_type, LambdaMV, gammaMV, q2min, q2max, Qs2min, Qs2max" << endl;
        return 1;
    }
    string gdist_type = argv[1];
    double LambdaMV = strtod(argv[2], NULL);
    double gammaMV = strtod(argv[3], NULL);
    double q2min = strtod(argv[4], NULL);
    double q2max = strtod(argv[5], NULL);
    double Qs2min = strtod(argv[6], NULL);
    double Qs2max = strtod(argv[7], NULL);
    if (gdist_type == "MV") {
        MVGluonDistribution gdist(LambdaMV, gammaMV, q2min, q2max, Qs2min, Qs2max);
        handle_input(gdist);
    }
    else if (gdist_type == "fMV") {
        FixedSaturationMVGluonDistribution gdist(LambdaMV, gammaMV, q2min, q2max, Qs2min);
        handle_input(gdist);
    }
    else {
        cerr << "Unrecognized gluon distribution type " << gdist_type << endl;
        return 1;
    }
    return 0;
}
#endif
