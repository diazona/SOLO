// Program for solving Balitsky-Kovchegov
// nonlinear equation 
// Options includ the running coupling in the Balitsky scenario and the parent dipole prescription
// Last changes 11.12.2012 AMS
// 

#include<cmath>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<ctime>
#include<cstdlib>
#include<cstring>
#include "comm_constr.h"
#include "input_funr.h"
#include "kovr.h"

using namespace std;

//**************************************************************************
double *h_new;	// vector that holds the result of iteration at Y+dY
double *h0;     // vector that holds the result of iteration at Y
double *h_old;  // vector that holds the result of itertation at Y+dY (need a second one)
const double xinit =  0.01; // initial value of x from which the evolution is started. This is the value that
// GBW model is evaluated and set in as the initial condition

const double rmin = exp(lrmin);		// the minimum value of r
const double rmax = exp(lrmax);		// the maximum value of r

double *xwezly;			// these are vectors needed for the Gaussian integration
double *wwezly;

double rG,logrG,lrsav;   // variables that store r and log(r) as global variables
double nxy0,nxy1;		// values of the amplitude evaluated at r

int QGWEZLY;			// parameter for Gaussian integration

// switch parameters, explained below
enum evtype {BK,BFKL};
enum RegType {CUTOFF,FROZEN};
enum startstatus {NEW,OLD};
enum alphtype {FIX,RUNBAL,RUNPAR,RUNMIN};
enum inptype { GBW, MODEL1 };
inptype inpmodel;
alphtype alphas;
evtype eqtype;
RegType Reg;
startstatus startup;

main()
{
    const int FieldWidth = 18;
    const int FieldWidthQ2 = 18;
    const int PrecisionNumber = 15;
    //
    // If we want to shift the initial rapidity by some constant
    // It does not affect the evolution
    const double Rapinitial = 0.0;
    
    clock_t  start,end;		// variables to measure time
    int i,j,k,iy;
    int CheckAccuracy=0;
    double logr;
    const int n=NR;			// size of the grid in r
    fstream fileout,fileout2;
    
    StartProgram();
    
    fileout2.open("kovr_fine_grid.dat",ios::out);
    
    // alpha strong
    // FIX fixed
    // RUNBAL running coupling Balitsky prescription
    // RUNPAR running coupling parent dipole prescription
    // RUNMIN running coupling minimum dipole prescription
    alphas = FIX;
    // linear or nonlinear equation
    // BFKL linear evolution
    // BK nonlinear evolution
    eqtype = BK;
    //what kind of initial condition
    // GBW 
    // MODEL1 some sample model
    inpmodel  = GBW;
    
    // initialisation of the vectors needed for the Gaussian integration
    
    double * xtemp = new double[NWEZLY*2];
    double * wtemp = new double[NWEZLY*2];
    xwezly = new double[NWEZLY];
    wwezly = new double[NWEZLY];
    
    gauleg(-1.0,1.0,xtemp,wtemp,2*NWEZLY);
    
    for(i=0;i<NWEZLY;i++)
    {
        xwezly[i] = xtemp[i+NWEZLY+1];
        wwezly[i] = wtemp[i+NWEZLY+1];
    }
    QGWEZLY = NWEZLY;
    
    // Booking the matrices
    Stars();
    cout << "Booking the matrices ..." << endl;
    
    h_new = new double[n];
    h_old = new double[n];
    h0 = new double[n];
    
    cout << "Done..." << endl;
    Stars();
    
    //Initialise the vectors with the input
    
    InitialiseMatrix(h0,NR);
    InitialiseMatrix(h_old,NR);
    
    // Print the vector to check 
    PrintMatrix(h0,NR);
    
    //  Zero the vector which will hold the result
    ZeroMatrix(h_new,NR);
    
    // Print some parameters just to check	
    PrintParameters();
    Stars();
    // cout << "Rapidity points : " << NY << endl;
    
    // Open the file to write
    fileout.open("kovr_coarse_grid.dat",ios::out);
    
    // Print out the inital values at rapidity=0
    // h0[i] and input(logr,xinit) should actually be the same
    
    for(i=0;i<NR;i++) {
        
        logr = lrmin + (double)i * Deltar;
        
        
        cout <<  Rapinitial << " " << exp(logr) << "  "  
        <<  input(logr,xinit) << " " << h0[i]  << endl;
        
        fileout <<  Rapinitial << " " << exp(logr) << "  "  
        <<  input(logr,xinit) << " " << h0[i]  << endl;
        
    }
    
    // write out the fine grid at Y=0
    {
        double res=0.0,resinput=0.0;
        for(int k=0;k<NN;k++) {
            logr = log(bmin + double(k)/double(NN)*(bmax-bmin));
            res = 1.0 - inter(logr,h_old);
            // double yout = 10.0;
            resinput = 1.0- input(logr,xinit);
            fileout2  << Rapinitial << " " << exp(logr) << "  " << setprecision(PrecisionNumber) << resinput << " "  <<  setprecision(PrecisionNumber) <<res <<  endl;
            
        }
    }
    
    // Loop over the rapidity points
    for(iy=0;iy<NY;iy++) {
        start = clock();
        
        // Loop over the iterations
        for(int iter=0;iter<IterMax;iter++)
        {
            if(iter==0&&iy==0) 
            {
                startup=NEW;
            }
            else
            {startup=OLD; }
            
            // call the function to solve this one step and perform the integrations
            SolveOneStep(NR);
            
            // copy vector h_new into h_old
            CopyMatrix(h_new,h_old,n);
            
        }
        //      if(CheckAccuracy==1) WriteOutErrors(h_new,h_old,NR);
        
        // write out the results for this rapidity to the file and onto screen
        double yout = Rapinitial+Deltay*(double)(iy+1);
        
        for(i=0;i<NR;i++) {
            
            logr = lrmin + (double)i * Deltar;
            
            
            cout <<  yout << " " << exp(logr) << "  "  
            <<  input(logr,xinit/exp(yout)) << " " << h_new[i]  << endl;
            
            fileout <<  yout << " " << exp(logr) << "  "  
            <<  input(logr,xinit/exp(yout)) << " " << h_new[i]  << endl;
        }
        // copy the matrix h_new into h0 before starting to move to another point in rapidity
        CopyMatrix(h_new,h0,NR);
        end = clock();
        cout << "Rapidity: " << yout << "    Time: " << (end-start)/CLOCKS_PER_SEC << " seconds "<< endl;
        
        if ((int(yout)-yout)==0) {
            // for integer values of rapidity, interpolate in r
            double res=0.0,resinput=0.0;
            for(int k=0;k<NN;k++)
            {
                logr = log(bmin + double(k)/double(NN)*(bmax-bmin));
                res = inter(logr,h_new);
                res = 1.0 - res;
                // double yout = 10.0;
                resinput = 1.0- input(logr,xinit/exp(yout));
                fileout2  << yout << " " << exp(logr) << "  " << setprecision(PrecisionNumber) << resinput << " "  <<  setprecision(PrecisionNumber) <<res <<  endl;
                
            }}
            
    }
    fileout.close();
    
    fileout2.close();
    Stars();
    cout << "Finishing the program. Freeing the memory ..." << endl;
    
    
    delete [] h_new ;
    delete [] h_old;
    delete [] h0 ;
    Stars();
    return 0;
}
//***************************************************************
// The input distribution, the initial condition for the Balitsky-Kovchegov equation

double input(double logr,double x)
{
    double Qs2;
    double r = exp(logr);
    
    switch (inpmodel) {
        // standard GBW model
        case GBW:
            return 1.0 - exp( - r * r * Qsat2(x) / 4.0);
            break;
            // some other model for testing
        case MODEL1:
            Qs2 = 0.01; 
            return r * r * Qs2 / 4.0 / ( 1.0 + r * r * Qs2 / 4.0);
            break;
    }	
}
//***************************************************************
double Qsat2(double x)
{
    // GBW model parametrization of the saturation scale	
    double lambda =  0.288;
    double Q02 = 0.56*pow(208, 1./3.);
    double x0 = 0.000304;
    
    return pow(x0/x,lambda) * Q02;
}
//***************************************************************
// function to initialise the vector with values from the input model
void InitialiseMatrix(double *h1,int n)
{
    int i;
    double logr;
    for(i=0;i<n;i++) {
        
        logr = lrmin + (double)i  * Deltar;
        h1[i] = input(logr,xinit);
    }	
}
//***************************************************************
void Stars()
{
    cout << "**********************************************" <<endl;
}
//***************************************************************
void StartProgram()
{
    Stars();
    cout << "*    Program solving                         *" << endl;
    cout << "*  Balitsky - Kovchegov equation             *" << endl;
    cout << "*           December 2012                    *" << endl;
    Stars();                                      
}
//***************************************************************
void PrintParameters()
{
    cout << " Parameters for this run: " << endl;
    cout << " r_min = " << rmin << "  r_max = " << rmax  << endl;
    cout << " Number of divisions: in r = " << NR << endl;
    cout << " Grid logarithmic in r " << endl;
    cout << " delta r = " << Deltar << endl;
    cout << " delta Y = " << Deltay << endl;
    cout << " Y_max   = " << Ymax << endl;
    cout << " asb      = " << alphastrong << endl;
    cout << " Number of iterations: " << IterMax << endl;
}
//***************************************************************
// Puts all vector elements to zero
void ZeroMatrix(double *h1,int n)
{
    for(int i=0;i<n;i++)
    {
        h1[i] = 0.0;
    }
}
//***************************************************************
// Copy matrix h1 into matrix h2
void CopyMatrix(double *h1,double *h2,int n)
{
    for(int i=0;i<n;i++){
        h2[i] = h1[i];
    }
}
//***************************************************************
// Print out matrix h1 onto screen
void PrintMatrix(double *h1,int n)
{
    for(int i=0;i<n;i++){
        cout << i <<  " " << exp(lrmin + (double)i * Deltar) << " " << h1[i] << " " << endl;
    }
}
//***************************************************************
void WriteOutErrors(double *h1,double *h2,int n)
{
    for(int i=0;i<n;i++){
        cout << i << 
        h1[i] << " " << h2[i] << "  " << 
        dabs(100.0*(h2[i]-h1[i])/(h2[i])) << endl;;
    }
}
//***************************************************************
// call the integration routine
void SolveOneStep(int n)
{
    int m = 5;  // this parameter controls the number of divisions of the integral in r
    // can be used to test accuracy of the calculation
    
    double deltam=(lrmax-lrmin)/(double)m;
    double sum = 0.0;
    
    // loop over all points in r grid
    for(int i=0;i<n;i++)
    {	
        logrG = lrmin + (double)i * Deltar;
        rG = exp(logrG);		
        nxy1 = h_old[i];
        nxy0 = h0[i];		
        sum = 0.0;
        
        for(int k=0;k<m;k++)
        {
            sum += quad2d(Kernel,lrmin+(double)k*deltam,lrmin+(double)(k+1)*deltam);
        }
        h_new[i] = h0[i] +  KERCOFF * sum;
    }
}
//***************************************************************
double qgauss(double (*func)(double),double a,double b)
{
    int j;
    double xr,xm,dx,s;
    
    static double x[]={0.076526521133497333755,0.227785851141645078080,
    0.373706088715419560673,0.510867001950827098004,0.636053680726515025453,
    0.746331906460150792614,0.839116971822218823395,0.912234428251325905868,
    0.963971927277913791268,0.993128599185094924786} ;
    static double w[]={0.152753387130725850698,0.149172986472603746788,
    0.142096109318382051329,0.131688638449176628898,0.118194531961518417312,
    0.101930119817240435037,0.083276741576704748725,0.062672048334109063570,
    0.040601429800386941331,0.017614007139152118312} ;
    xm=0.5*(a+b);
    xr=0.5*(b-a);
    s=0;
    for(j=0;j<10;j++)
    {
        dx=xr*x[j];
        s+=w[j]*((*func)(xm+dx)+(*func)(xm-dx));
    }
    return s*=xr;
}
// *********************************************************************
double qgauss24(double (*func)(double),double a,double b)
{
    int j;
    double xr,xm,dx,s;
    static double x[]={0.032380,0.097005,0.161222,0.224764,0.287362,
    0.348756,0.408686,0.466903,0.523161,0.577225,0.628867,0.677872,
    0.724034,0.767159,0.807066,0.843588,0.876572,
    0.905879,0.931387,0.952988,0.970592,0.984125,0.993530,0.998771} ;
    static double w[]={0.064738,0.064466,0.063924,0.063114,0.062039,
    0.060704,0.059115,0.057277,0.055200,0.052890,0.050359,0.047617,
    0.044675,0.041545,0.038241,0.034777,0.031167,0.027427,0.023571,
    0.019616,0.015579,0.011477,0.007328,0.003153};
    
    xm=0.5*(a+b);
    xr=0.5*(b-a);
    s=0;
    for(j=0;j<24;j++)
    {
        dx=xr*x[j];
        s+=w[j]*((*func)(xm+dx)+(*func)(xm-dx));
    }
    return s*=xr;
}
//***************************************************************
/* this routine calculates double integral */
/* together with the qgauss, w1,w2 functions */
double quad2d(double (*func)(double,double),double lrlow,double lrhi)
{
    ur2func=func;	
    return QGAUSS(w1,lrlow,lrhi);
}
//***************************************************************
/* this is the interface to the double integration */
// in this function we actually compute the integral over the angle
double w1(double lrprim)
{
    int m = 1;
    lrsav = lrprim;
    double deltam = 2.0* M_PI /(double)m;
    double sum = 0.0;
    
    for(int k=0;k<m;k++)
    {
        sum+= QGAUSS(w2,(double)k*deltam,(double)(k+1)*deltam);
    }
    return sum;
}
// **************************************************************
double w2(double angle)
{
    return (*ur2func)(lrsav,angle);	
}
//***************************************************************
// Main routine with the kernel
// Last changes 14.11.2012 A.S.
double Kernel(double lrprim,double phi)
{ 
    double kernelint;
    double nzx0,nzx1,nzy0,nzy1;
    enum stat_bound {YES,NO};
    stat_bound SetZero=NO;
    double eps = 0.00001;
    // parsprim contains information about dipole XZ 
    double parsprim;
    // parsbis contains information about dipole YZ 
    double parsbis;
    
    double rprim = exp(lrprim);
    double lrbis ; 
    
    parsprim = lrprim;
    
    //Second dipole ZY
    //Size of ZY
    double rbis = sqrt(rG * rG + rprim * rprim + 2.0 * rG * rprim * cos(phi)) ;
    
    // is this value inside the grid
    if(rbis<rmin) {rbis=rmin;SetZero=YES;}
    else if(rbis>rmax) {rbis=rmax;SetZero=YES;}		
    lrbis = log(rbis);
    parsbis=lrbis;
    
    // if not the first step interpolate the vectors
    if(startup==OLD)
    { 
        nzx0  = inter(parsprim,h0);
        nzy0  = inter(parsbis,h0);
        nzx1  = inter(parsprim,h_old);
        nzy1  = inter(parsbis,h_old);
    }
    
    else 
    {
        // if the first step call directly the input
        nzx0  = input(parsprim,xinit);
        nzy0  = input(parsbis,xinit);
        nzx1  = input(parsprim,xinit);
        nzy1  = input(parsbis,xinit);		
    }
    // kernelint is alphas * K * rprim * rprim
    // K is the BFKL kernel
    // rprim *rprim comes from jacobian and the fact that we integrated over log(r)
    
    switch (alphas) {
        // fixed coupling constant
        case FIX:
            kernelint = alphastrong * rG * rG  / (rbis * rbis);
            break;
            // running coupling using Balitsky prescription
        case RUNBAL:
            kernelint = alphalo(rG) * ((alphalo(rprim)/alphalo(rbis)-1.0) +
            (rprim*rprim)/(rbis*rbis)*(alphalo(rbis)/alphalo(rprim)-1.0)+rG*rG/(rbis*rbis) );
            break;
            // running coupling with argument of a parent dipole
        case RUNPAR:
            kernelint = alphalo(rG) * rG * rG  / (rbis * rbis);
            break;
            // running coupling with argument of a minimum size dipole
        case RUNMIN:
            kernelint = alphalo(dmin(rG,dmin(rprim,rbis))) * rG * rG  / (rbis * rbis);
            break;	
    }
    
    
    if(Reg==CUTOFF&&SetZero==YES) {
        return 0.0;
    }
    else{
        switch(eqtype){
            // This is for solution of nonlinear equation
            case BK:
                return  kernelint * 
                ( 0.5 * ( nzx0 + nzx1 + nzy0 + nzy1 - nxy0 - nxy1 )  -
                1.0/3.0 * ( nzx0 * nzy0 + nzx1 * nzy1 ) -
                1.0/6.0 * ( nzx0 * nzy1 + nzx1 * nzy0 ));
                break;
                // This is for solution of linear equation
            case BFKL:
                return  kernelint *
                0.5 * ( nzx0 + nzx1 + nzy0 + nzy1 - nxy0 - nxy1 );
                break;
        }
    }
}
// **************************************************************
double QGAUSS(double (*func)(double),double a,double b)
{
    double dx;
    
    double xm=0.5*(a+b);
    double xr=0.5*(b-a);
    double s=0;
    for(int j=0;j<QGWEZLY;j++)
    {
        dx=xr*xwezly[j];
        s+=wwezly[j]*((*func)(xm+dx)+(*func)(xm-dx));
    }
    return s*=xr;
}
// **************************************************************
double dabs(double arg)
{
    if(arg>=0.0) return arg;
                                                else return -arg;
}
//***************************************************************
double dmax(double arg1,double arg2)
{
    if(arg1>=arg2) return arg1;
                                                else return arg2;
}
//***************************************************************
double dmin(double arg1,double arg2)
{
    if(arg1>=arg2) return arg2;
                                                else return arg1;
}
//***************************************************************
// alpha strong in LO. contains freezing at scale mu2
double alphalo(double arg)
{
    double beta = (33.0 - 2.0 * NF) / (12.0 * M_PI);
    return NC / M_PI / beta/ log((1.0/(arg * arg) + mu2) / Lambda2)	;
}
//***************************************************************



