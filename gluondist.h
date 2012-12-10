#ifndef _GLUONDIST_H_
#define _GLUONDIST_H_

class GluonDistribution {
public:
    virtual double S2(double r2, double Qs2) = 0;
    virtual double S4(double r2, double s2, double t2, double Qs2) = 0;
};

class GBWGluonDistribution: public GluonDistribution {
public:
    double S2(double r2, double Qs2);
    double S4(double r2, double s2, double t2, double Qs2);
};

class MVGluonDistribution: public GluonDistribution {
public:
    MVGluonDistribution(double LambdaMV) : LambdaMV(LambdaMV) {};
    double S2(double r2, double Qs2);
    double S4(double r2, double s2, double t2, double Qs2);
private:
    double LambdaMV;
};

#endif // _GLUONDIST_H_