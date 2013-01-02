#include <cassert>
#include <iostream>
#include <cmath>
#include "dss_pinlo.h"

using namespace std;

int main(int argc, char** argv) {
    DSSpiNLO ffs("PINLO.DAT");
    double zlist[] = {exp(-4.605), exp(-3.912), 0.8};
    double qs2list[] = {exp(0), exp(0.2231), 3.91213};
    for (int iz = 0; iz < sizeof(zlist) / sizeof(zlist[0]); iz++) {
        for (int iq = 0; iq < sizeof(qs2list) / sizeof(qs2list[0]); iq++) {
            double z = zlist[iz];
            double qs2 = qs2list[iq];
            ffs.update(z, qs2);
            cout << "f_up/p+(" << z << "," << qs2 << ") = " << ffs.fragmentation(DSSpiNLO::up, DSSpiNLO::pi_plus) << endl;
            cout << "f_ubar/p+(" << z << "," << qs2 << ") = " << ffs.fragmentation(DSSpiNLO::up_bar, DSSpiNLO::pi_plus) << endl;
            cout << "f_down/p+(" << z << "," << qs2 << ") = " << ffs.fragmentation(DSSpiNLO::down, DSSpiNLO::pi_plus) << endl;
            cout << "f_dbar/p+(" << z << "," << qs2 << ") = " << ffs.fragmentation(DSSpiNLO::down_bar, DSSpiNLO::pi_plus) << endl;
            cout << "f_g/p+(" << z << "," << qs2 << ") = " << ffs.fragmentation(DSSpiNLO::gluon, DSSpiNLO::pi_plus) << endl;
        }
    }
    return 0;
}