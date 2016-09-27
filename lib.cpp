#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>

#include <gsl/gsl_integration.h>

using namespace std;

#include "lib.h"

double integrate_1D(double (*f)(double), double * xpts, double * xwts, int nx)
{
    double sum = 0.0;
    for (int ix = 0; ix < nx; ix++)
        sum += xwts[ix] * (*f)(xpts[ix]);

    return (sum);
}

complex<double> integrate_1D_FT(double (*f)(double), double * xpts, double * xwts, int nx, double k)
{
    complex<double> sum = (0,0);
	for (int ix = 0; ix < nx; ix++)
		sum += xwts[ix] * exp(- i * k * xpts[ix]) * (*f)(xpts[ix]);

    return (sum);
}


double integrate_2D(double (*f)(double, double), double * xpts, double * ypts, double * xwts, double * ywts, int nx, int ny)
{
	double sum = 0.0;
	for (int ix = 0; ix < nx; ix++)
	for (int iy = 0; iy < ny; iy++)
		sum += xwts[ix] * ywts[iy] * (*f)(xpts[ix], ypts[iy]);

	return (sum);
}

double integrate_1D(double (*f)(double), double Lx, int nx)
{
    double sum = 0.0;
    for (int ix = 1; ix < nx; ix++)
	{
		double xp = (double)ix * M_PI / (double)nx;
		double s1 = sin(xp);
        sum += (Lx * M_PI) * (*f)(Lx * cot(xp)) / (nx*s1*s1);
	}

    return (sum);
}


double integrate_2D(double (*f)(double, double), double Lx, double Ly, int nx, int ny)
{
    double sum = 0.0;
    for (int ix = 1; ix < nx; ix++)
    for (int iy = 1; iy < ny; iy++)
	{
		double xp = (double)ix * M_PI / (double)nx;
		double yp = (double)iy * M_PI / (double)ny;
		double s1 = sin(xp);
		double s2 = sin(yp);
        sum += (Lx * Ly * M_PI * M_PI) * (*f)(Lx * cot(xp), Ly * cot(yp)) / (nx*ny*s1*s1*s2*s2);
	}

    return (sum);
}

