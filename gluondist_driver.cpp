#include <cstdlib>
#include <iostream>
#include "context.h"
#include "gluondist.h"

// http://stackoverflow.com/a/1404473/56541
// use "extern" explicitly to achieve external linkage
extern const double inf = 10;
ostream& logger = cerr;

void momentum_space_query(GluonDistribution* gdist) {
    double q2, Y;
    cout << "q2\tY\tx\tQs2\tF" << endl;
    while (cin >> q2 >> Y) {
        cout << q2 << "\t"
                << Y << "\t"
                << exp(-Y) << "\t"
                << gdist->Qs2(Y) << "\t"
                << gdist->F(q2, Y) << endl;
    }
}

void position_space_query(GluonDistribution* gdist) {
    double r2, Y;
    cout << "r2\tY\tx\tQs2\tS" << endl;
    while (cin >> r2 >> Y) {
        cout << r2 << "\t"
                << Y << "\t"
                << exp(-Y) << "\t"
                << gdist->Qs2(Y) << "\t"
                << gdist->S2(r2, Y) << endl;
    }
}

void satscale_query(GluonDistribution* gdist) {
    double Y;
    cout << "Y\tx\tQs2" << endl;
    while (cin >> Y) {
        cout << Y << "\t"
             << exp(-Y) << "\t"
             << gdist->Qs2(Y) << endl;
    }
}

#define usage() cerr << "Usage: " << argv[0] << " <query|printgrid> <S2|F|Qs2> <filename.cfg>" << endl; return 1;


/**
 * A driver program that constructs a gluon distribution object
 * and does various things with it.
 */
int main(int argc, char** argv) {
    if (argc < 4) {
        usage();
    }
    string mode(argv[1]);
    string quantity(argv[2]);
    string filename(argv[3]);
    
    ContextCollection cc(filename);
    GluonDistribution* gdist;
    try {
        cc.create_contexts();
        gdist = cc[0].gdist;
    }
    catch (const exception& e) {
        cout << "Error in parsing: " << e.what() << endl;
        return 1;
    }
    
    if (mode == "printgrid") {
        try {
            if (quantity == "F") {
                gdist->write_pspace_grid(cout);
            }
            else if (quantity == "S2") {
                gdist->write_rspace_grid(cout);
            }
            else if (quantity == "Qs2") {
                gdist->write_satscale_grid(cout);
            }
        }
        catch (NoGridException& e) {
            cerr << "No grid available" << endl;
            return 1;
        }
        return 0;
    }
    else if (mode == "query") {
        if (quantity == "F") {
            momentum_space_query(gdist);
        }
        else if (quantity == "S2") {
            position_space_query(gdist);
        }
        else if (quantity == "Qs2") {
            satscale_query(gdist);
        }
    }
    return 0;
}