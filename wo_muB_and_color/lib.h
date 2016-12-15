#ifndef LIB_H
#define LIB_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>

#include <gsl/gsl_integration.h>

using namespace std;

const std::complex<double> i(0, 1);
inline double cot(double x){return (cos(x)/sin(x));}

double integrate_1D(double (*f)(double), double * xpts, double * xwts, int nx);
complex<double> integrate_1D_FT(double (*f)(double), double * xpts, double * xwts, int nx, double k);
double integrate_2D(double (*f)(double, double), double * xpts, double * ypts, double * xwts, double * ywts, int nx, int ny);
//double integrate_1D(double (*f)(double), double Lx, int nx);
//double integrate_2D(double (*f)(double, double), double Lx, double Ly, int nx, int ny);

//main interpolation routine
double interpolate1D(double * x, double * y, double x0, long size, int kind, bool uniform_spacing,
			bool returnflag = false, double default_return_value = 0.0);

//subsidiary interpolation routines
double interpLinearDirect(double * x, double * y, double x0, long size,
				bool returnflag = false, double default_return_value = 0.0);
double interpLinearNondirect(double * x, double * y, double x0, long size,
				bool returnflag = false, double default_return_value = 0.0);
double interpCubicDirect(double * x, double * y, double x0, long size,
				bool returnflag = false, double default_return_value = 0.0);
double interpCubicNonDirect(double * x, double * y, double x0, long size,
				bool returnflag = false, double default_return_value = 0.0);
long binarySearch(double * A, int length, double value,
				bool skip_out_of_range = true, bool verbose = false);

#endif
