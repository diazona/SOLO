#include<cmath>
#include<iostream>
#include<cstdlib>
#include "comm_constr.h"
#include "input_funr.h"

using namespace std;

//********************************************************
// Linear interpolation procedure in 1 dimension
//
double inter(double params,const double *h)
{

    int i,j,k;
    double t,u,w,t1,u1,w1;
    static double eps=1.0e-8;
    if(params>=lrmax)
    {
        //cout << "Inter: r " << exp(params[0]) << endl;
        params=lrmax-eps;
    }

    if(params<=lrmin)
    {
        //cout << "Inter: r " << exp(params[0]) << endl;
        params=lrmin+eps;
    }

    i = (int)((params-lrmin)/Deltar);

    if(i<0||i>=NR)
    {   cout << i  << endl;
        exit(0);
    }
    t = (params-lrmin)/Deltar - (double)i;

    t1 = 1.0 - t;


    return t1 *  h[i]
           +  t *  h[i+1] ;

}
//*****************************************************************
