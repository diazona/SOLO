#include <cmath>
#include "gluondist.h"

double GBWGluonDistribution::S2(double r2, double Qs2) {
    return exp(-0.25 * r2 * Qs2);
}
double GBWGluonDistribution::S4(double r2, double s2, double t2, double Qs2) {
    return exp(-0.25 * Qs2 * (s2 + t2));
}

double MVGluonDistribution::S2(double r2, double Qs2) {
    return pow(M_E + 1.0 / (sqrt(r2) * LambdaMV), -0.25 * r2 * Qs2);
}
double MVGluonDistribution::S4(double r2, double s2, double t2, double Qs2) {
    return S2(s2, Qs2) * S2(t2, Qs2);
}
