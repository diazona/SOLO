#ifndef _COUPLING_H_
#define _COUPLING_H_

#include <cmath>

class Coupling {
public:
    virtual double alphasbar(double kT2) = 0;
};

class FixedCoupling : public Coupling {
private:
    double value;
public:
    FixedCoupling(double alphasbar) : value(alphasbar) {}
    double alphasbar(double kT2) {
        return value;
    }
};

class LORunningCoupling : public Coupling {
private:
    double log_LambdaQCD;
    double inverse_beta_2pi;
    double regulator; // position of the Landau pole
public:
    LORunningCoupling(double LambdaQCD, double beta, double regulator) : log_LambdaQCD(log(LambdaQCD)), inverse_beta_2pi(0.5 / (M_PI * beta)), regulator(regulator) {}
    double alphasbar(double kT2) {
        return inverse_beta_2pi / (log(kT2 + regulator) - log_LambdaQCD);
    }
};

#endif // _COUPLING_H_