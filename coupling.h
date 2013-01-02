#ifndef _COUPLING_H_
#define _COUPLING_H_

#include <cmath>
#include <ostream>
#include <string>

class Coupling {
public:
    virtual double alphasbar(double kT2) = 0;
    virtual const char* name() = 0;
};

class FixedCoupling : public Coupling {
private:
    double value;
    std::string _name;
public:
    FixedCoupling(double alphasbar);
    double alphasbar(double kT2);
    const char* name();
};

class LORunningCoupling : public Coupling {
private:
    double log_LambdaQCD;
    double inverse_beta_2pi;
    double regulator; // position of the Landau pole
    std::string _name;
public:
    LORunningCoupling(double LambdaQCD, double beta, double regulator);
    double alphasbar(double kT2);
    const char* name();
};

std::ostream& operator<<(std::ostream& out, Coupling& cpl);

#endif // _COUPLING_H_
