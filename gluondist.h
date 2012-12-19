#ifndef _GLUONDIST_H_
#define _GLUONDIST_H_

#include "interp2d.h"

class GluonDistribution {
public:
    virtual double S2(double r2, double Qs2) = 0;
    virtual double S4(double r2, double s2, double t2, double Qs2) = 0;
    virtual double F(double q2, double Qs2) = 0;
};

class GBWGluonDistribution: public GluonDistribution {
public:
    double S2(double r2, double Qs2);
    double S4(double r2, double s2, double t2, double Qs2);
    double F(double q2, double Qs2);
};

class MVGluonDistribution: public GluonDistribution {
public:
    MVGluonDistribution(double LambdaMV, double q2min, double q2max, double Qs2min, double Qs2max);
    ~MVGluonDistribution();
    double S2(double r2, double Qs2);
    double S4(double r2, double s2, double t2, double Qs2);
    double F(double q2, double Qs2);
private:
    double LambdaMV;
    double q2min, q2max;
    double Qs2min, Qs2max;
    double* log_q_values;
    double* log_Qs2_values;
    double* F_dist;
    double* F_error;
    interp2d* interp_dist;
    interp2d* interp_error;
    gsl_interp_accel* q_accel;
    gsl_interp_accel* Qs2_accel;
};

#endif // _GLUONDIST_H_