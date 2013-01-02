#include <sstream>
#include "coupling.h"

using namespace std;

FixedCoupling::FixedCoupling(double alphasbar) : value(alphasbar) {
    ostringstream s;
    s << "Fixed(alphasbar = " << value << ")";
    _name = s.str();
}
double FixedCoupling::alphasbar(double kT2) {
    return value;
}
const char* FixedCoupling::name() {
    return _name.c_str();
}

LORunningCoupling::LORunningCoupling(double LambdaQCD, double beta, double regulator) :
    log_LambdaQCD(log(LambdaQCD)),
    inverse_beta_2pi(0.5 / (M_PI * beta)),
    regulator(regulator) {
    ostringstream s;
    s << "LORunning(LambdaQCD = " << LambdaQCD << ", beta = " << beta << ", regulator = " << regulator << ")";
    _name = s.str();
}
double LORunningCoupling::alphasbar(double kT2) {
    return inverse_beta_2pi / (log(kT2 + regulator) - log_LambdaQCD);
}
const char* LORunningCoupling::name() {
    return _name.c_str();
}


std::ostream& operator<<(std::ostream& out, Coupling& cpl) {
    out << cpl.name();
    return out;
}
