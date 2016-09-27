#ifndef DEFS_H
#define DEFS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>

//#include <gsl/gsl_integration.h>

using namespace std;

#include "lib.h"

/*USAGE: debugger(__LINE__, __FILE__);*/
void inline debugger(int cln, const char* cfn)
{
	cerr << "You made it to " << cfn << ":" << cln << "!" << endl;
	return;
}

extern double vs, Neff, tauf, tau0, Tf, T0, nu, nuVB, ds, A, m;
extern double mByT, alpha, alpha0, phi0;

extern const int n_xi_pts;
extern const int n_k_pts;
extern const int n_tau_pts;
extern double * xi_pts_0_inf, * xi_wts_0_inf, * xi_pts_minf_inf, * xi_wts_minf_inf;
extern double * k_pts, * k_wts;
extern double * tau_pts, * tau_wts;

// my additional functions
/*double zeta(double x1, double x2);
double gamma(double x1, double x2);
double krho(double x);
double komega(double x);
double alpha_int(double x1, double x2);
double phi0_int(double x1, double x2);
double sym(double x, double xp);
double phi(double y, double x, double xp);
double phip(double y, double x, double xp);
double phim(double y, double x, double xp);
double phipp(double y, double x, double xp);
double phimp(double y, double x, double xp);
double phipp00(double xi1, double xi2);
double phipp20(double xi1, double xi2);
double phipp02(double xi1, double xi2);
double phimp11(double xi1, double xi2);*/

inline double Temperature(double tau)
{
	return (T0 * pow(tau0 / tau, vs*vs));
}

inline double dTemperaturedtau(double tau)
{
	return (- T0 * vs * vs * pow(tau0 / tau, vs*vs) / tau);
}

inline double sEnt(double tau)
{
	return (2.0 *M_PI * M_PI *Neff/45.0)*pow(Temperature(tau), 3.0);

}

inline double wEnth(double tau)
{
	return (sEnt(tau) * Temperature(tau));
}

inline double invTemperature(double Tf)
{
	return ( tauf*pow( Tf/T0, 1./(vs*vs) ) );
}

inline complex<double> beta(double k)
{
	return (std::sqrt(std::complex<double>(alpha*alpha - vs*vs*k*k)));
}

inline double incompleteGamma3(double x)
{
	return (exp(-x)*(2.+x*(2.+x)));
}

inline double incompleteGamma4(double x)
{
	return (exp(-x)*(6. + x*(6. + x*( 3. + x))));
}

inline double incompleteGamma5(double x)
{
	return (exp(-x)*(24.+x*(24.+x*(12.+x*(4.+x)))));
}

inline double norm_int(double x)
{
	double cx = cosh(x);
	return (incompleteGamma3(mByT * cx) / (cx*cx));
}

inline double Frho(double x)
{
	double cx = cosh(x);
	return (vs*vs*incompleteGamma4(mByT * cx) / (cx*cx));
}

inline double Fomega(double x)
{
	double cx = cosh(x);
	return (tanh(x)*incompleteGamma4(mByT * cx) / (cx*cx));
}

inline complex<double> Ftilde_rho(double k)
{
	return (integrate_1D_FT(Frho, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts, k));
}

inline complex<double> Ftilde_omega(double k)
{
	return (integrate_1D_FT(Fomega, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts, k));
}

inline complex<double> Gtilde_rho(double k, double t, double t_p)
{
	double t_by_tp = t / t_p;
	double t_x_T_at_t = t * Temperature(t);
	double tp_x_T_at_tp = t_p * Temperature(t_p);
	double log_tBYtp = log(t_by_tp);
	double prefactor = pow(t_by_tp, -alpha);
	complex<double> b = beta(k);
	complex<double> factor1 = cosh(b * log_tBYtp) + ((alpha + k*k + (nuVB*k*k / (2.*tp_x_T_at_tp))) / b) * sinh(b * log_tBYtp);
	double factor2 = exp(-(nuVB*k*k / (4.0 * alpha)) *((1.0/tp_x_T_at_tp) - (1.0/t_x_T_at_t)));

	return (prefactor * factor1 * factor2);
}

inline complex<double> Gtilde_omega(double k, double t, double t_p)
{
	double t_by_tp = t / t_p;
	double T_at_t = Temperature(t);
	double t_x_T_at_t = t * T_at_t;
	double T_at_tp = Temperature(t_p);
	double tp_x_T_at_tp = t_p * T_at_tp;
	double log_tBYtp = log(t_by_tp);

	double a = alpha;
	complex<double> b = beta(k);

	complex<double> prefactor = i * pow(t_by_tp, 1.0-a) / (8.0 * k * a * b * t_x_T_at_t * t_x_T_at_t * T_at_tp);
	complex<double> factor1a = 2.0*k*k*b * cosh(b * log_tBYtp) * (-nuVB * T_at_t * tp_x_T_at_tp
																	+ 2.0 * a * t * T_at_t * T_at_t * (nuVB + 2.0 * tp_x_T_at_tp)
																	- nuVB * t * tp_x_T_at_tp * dTemperaturedtau(t));
	complex<double> factor1b = sinh(b * log_tBYtp) * (
								T_at_t * ( ( 8.0*a*t_x_T_at_t*tp_x_T_at_tp*b*b ) - ( k*k*nuVB+4.0*a*a*t_x_T_at_t ) * ( k*k*nuVB+2.0*(k*k + a)*tp_x_T_at_tp )  )
								- k*k*nuVB*t*( k*k*nuVB + 2.0*(k*k + a)*tp_x_T_at_tp )*dTemperaturedtau(t)
								);
	complex<double> factor1 = factor1a + factor1b;
	double factor2 = exp(-(nuVB*k*k / (4.0 * a)) *((1.0/tp_x_T_at_tp) - (1.0/t_x_T_at_t)));
	
	return (prefactor * factor1 * factor2);
}

inline complex<double> Ctilde_rho_rho(double k)
{
	double pref = 2.0 * nu / A;
    complex<double> sum = 0.0;
    for (int it = 0; it < n_tau_pts; it++)
	{
		double tau = tau_pts[it];
		complex<double> Gtr = Gtilde_rho(k, tauf, tau);
        sum += tau_wts[it] * Gtr * conj(Gtr) / (tau*tau*tau*wEnth(tau));
	}

    return (pref * sum);
}

inline complex<double> Ctilde_rho_omega(double k)
{
	double pref = 2.0 * nu / A;
    complex<double> sum = 0.0;
    for (int it = 0; it < n_tau_pts; it++)
	{
		double tau = tau_pts[it];
		complex<double> Gtr = Gtilde_rho(k, tauf, tau);
		complex<double> cGto = Gtilde_omega(-k, tauf, tau);
        sum += tau_wts[it] * Gtr * cGto / (tau*tau*tau*wEnth(tau));
	}

    return (pref * sum);
}

inline complex<double> Ctilde_omega_rho(double k)
{
	double pref = 2.0 * nu / A;
    complex<double> sum = 0.0;
    for (int it = 0; it < n_tau_pts; it++)
	{
		double tau = tau_pts[it];
		complex<double> Gto = Gtilde_omega(k, tauf, tau);
		complex<double> cGtr = Gtilde_rho(-k, tauf, tau);
        sum += tau_wts[it] * Gto * cGtr / (tau*tau*tau*wEnth(tau));
	}

    return (pref * sum);
}

inline complex<double> Ctilde_omega_omega(double k)
{
	double pref = 2.0 * nu / A;
    complex<double> sum = 0.0;
    for (int it = 0; it < n_tau_pts; it++)
	{
		double tau = tau_pts[it];
		complex<double> Gto = Gtilde_omega(k, tauf, tau);
        sum += tau_wts[it] * Gto * conj(Gto) / (tau*tau*tau*wEnth(tau));
	}

    return (pref * sum);
}

////////////////////////////////////////////////////////////////////////////////
// All functions below this point correspond to HBT correlator
////////////////////////////////////////////////////////////////////////////////

inline double zeta(double x1, double x2)
{
    return ( Tf / (cosh(x1) + cosh(x2)) );
}

inline double gamma(double x1, double x2)
{
    double shx1_minus_shx2 = sinh(x1) - sinh(x2);
    return (shx1_minus_shx2 * shx1_minus_shx2);
}

inline double krho(double x)
{
    double chx = cosh(x);
    return ( (vs*vs / Tf) * chx * chx );
}

inline double komega(double x)
{
    return ( cosh(x) * sinh(x) );
}

inline double alpha_int(double x1, double x2)
{
    double z = zeta(x1,x2);
    double g = incompleteGamma4(m/z);
    return (cosh(x1) * cosh(x2) * pow(z,4.0) * g);
}

inline double phi0_int(double x1, double x2)
{
    double z = zeta(x1, x2);
    double g = gamma(x1,x2);
    double g4 = incompleteGamma4(m/z);
    return (-2.0 * alpha0 * g * cosh(x1) * cosh(x2) * pow(z, 4.0) * g4);
}

inline double get_phi0(double y)
{
	double sum = 0.0;
	for (int ix1 = 0; ix1 < n_xi_pts; ++ix1)
	for (int ix2 = 0; ix2 < n_xi_pts; ++ix2)
	{
		double x1 = xi_pts_minf_inf[ix1];
		double x2 = xi_pts_minf_inf[ix2];
		double z = zeta(x1, x2);
		double g = gamma(y-x1,y-x2);
		double g4 = incompleteGamma4(m/z);
		sum += -2.0 * xi_wts_minf_inf[ix1] * xi_wts_minf_inf[ix2] * alpha0 * g * cosh(x1) * cosh(x2) * pow(z, 4.0) * g4;
	}
	return (sum);
}

inline double sym(double x, double xp)
{
    double z_x_xp = zeta(x, xp);
    return (cosh(xp) * pow(z_x_xp, 5.0) * incompleteGamma5(m/z_x_xp) );
}

inline double phi(double y, double x, double xp)
{
    double g = gamma(y - x, y - xp);
    return (alpha0 * sym(x, xp) * (g + phi0) );
}

inline double phip(double y, double x, double xp)
{
    return (0.5 * (phi(y, x, xp) + phi(y, -x, xp)));
}

inline double phim(double y, double x, double xp)
{
    return (0.5 * (phi(y, x, xp) - phi(y, x, xp)));
}

inline double phipp(double y, double x, double xp)
{
    return (0.5 * (phip(y, x, xp) + phip(y, x, -xp)));
}

inline double phimp(double y, double x, double xp)
{
    return (0.5 * (phim(y, x, xp) + phim(y, x, -xp)));
}

inline double phipp00(double xi1, double xi2)
{
    return ( alpha0 * (-1.0 + phi0 + cosh(xi1) * cosh(xi2)) * sym(xi1, xi2) );
}

inline double phipp20(double xi1, double xi2)
{
    return (0.5 * alpha0 * (cosh(2.* xi1) - 2.* cosh(xi1) * cosh(xi2) + cosh(2. * xi2)) * sym(xi1, xi2) );
}

inline double phipp02(double xi1, double xi2)
{
    return (0.5 * alpha0 * (cosh(2. * xi1) - 2.* cosh(xi1) * cosh(xi2) + cosh(2.* xi2)) * sym(xi1, xi2) );
}

inline double phimp11(double xi1, double xi2)
{
    return (2.* alpha0 * (-cosh(xi1) + cosh(xi2)) * sinh(xi1) * sym(xi1, xi2) );
}

inline complex<double> Fbt_rho_pp00(double y, double k)
{
	double sum = 0.0;
	for (int ix = 0; ix < n_xi_pts; ix++)
	for (int iy = 0; iy < n_xi_pts; iy++)
	{
		double x = xi_pts_0_inf[ix];
		double y = xi_pts_0_inf[iy];
		sum += xi_wts_0_inf[ix] * xi_wts_0_inf[iy] * krho(x) * cos(k * x) * phipp00(x, y);
		/*double kr = krho(x);
		double pxy = phipp00(x, y);
		double ckx = cos(k * x);
		double sign = sgn(kr * ckx * pxy);
		double log_factor = log(abs(kr)+1e-100) + log(abs(ckx)+1e-100) + log(abs(pxy)+1e-100);
		sum += sign * xi_wts_0_inf[ix] * xi_wts_0_inf[iy] * exp(log_factor);*/
		//cout << "Fbt_rho_pp00(): " << x << "   " << y << "   " << krho(x) << "   " << cos(k * x) << "   " << phipp00(x, y) << endl;
	}
	return ( 4.0 * exp(-i * k * y) * sum );
}

inline complex<double> Fbt_rho_pp20(double y, double k)
{
	double sum = 0.0;
	for (int ix = 0; ix < n_xi_pts; ix++)
	for (int iy = 0; iy < n_xi_pts; iy++)
	{
		double x = xi_pts_0_inf[ix];
		double y = xi_pts_0_inf[iy];
		sum += xi_wts_0_inf[ix] * xi_wts_0_inf[iy] * krho(x) * cos(k * x) * phipp20(x, y);
		/*double kr = krho(x);
		double pxy = phipp00(x, y);
		double ckx = cos(k * x);
		double cy = cosh(y);
		double sign = sgn(kr * ckx * pxy);
		double log_factor = log(abs(kr)+1e-100) + log(abs(ckx)+1e-100) + log(abs(pxy)+1e-100) + 2.0*log(abs(cy)+1e-100);
		sum += sign * xi_wts_0_inf[ix] * xi_wts_0_inf[iy] * exp(log_factor);*/
		//cout << "Fbt_rho_pp20(): " << x << "   " << y << "   " << krho(x) << "   " << cos(k * x) << "   " << phipp20(x, y) << endl;
	}
	return ( 4.0 * exp(-i * k * y) * sum );
}

inline complex<double> Fbt_rho_pp02(double y, double k)
{
	double sum = 0.0;
	for (int ix = 0; ix < n_xi_pts; ix++)
	for (int iy = 0; iy < n_xi_pts; iy++)
	{
		double x = xi_pts_0_inf[ix];
		double y = xi_pts_0_inf[iy];
		sum += xi_wts_0_inf[ix] * xi_wts_0_inf[iy] * krho(x) * cos(k * x) * phipp02(x, y);
		/*double kr = krho(x);
		double pxy = phipp00(x, y);
		double ckx = cos(k * x);
		double sy = sinh(y);
		double sign = sgn(kr * ckx * pxy);
		double log_factor = log(abs(kr)+1e-100) + log(abs(ckx)+1e-100) + log(abs(pxy)+1e-100) + 2.0*log(abs(sy)+1e-100);
		sum += sign * xi_wts_0_inf[ix] * xi_wts_0_inf[iy] * exp(log_factor);*/
		//cout << "Fbt_rho_pp02(): " << x << "   " << y << "   " << krho(x) << "   " << cos(k * x) << "   " << phipp02(x, y) << endl;
	}
	return ( 4.0 * exp(-i * k * y) * sum );
}

inline complex<double> Fbt_rho_mp11(double y, double k)
{
	double sum = 0.0;
	for (int ix = 0; ix < n_xi_pts; ix++)
	for (int iy = 0; iy < n_xi_pts; iy++)
	{
		double x = xi_pts_0_inf[ix];
		double y = xi_pts_0_inf[iy];
		sum += xi_wts_0_inf[ix] * xi_wts_0_inf[iy] * krho(x) * sin(k * x) * phimp11(x, y);
		/*double kr = krho(x);
		double pxy = phipp00(x, y);
		double skx = sin(k * x);
		double sy = sinh(y);
		double cy = cosh(y);
		double sign = sgn(kr * skx * pxy * sy * cy);
		double log_factor = log(abs(kr)+1e-100) + log(abs(skx)+1e-100) + log(abs(pxy)+1e-100) + log(abs(cy)+1e-100);
		sum += sign * xi_wts_0_inf[ix] * xi_wts_0_inf[iy] * exp(log_factor);*/
		//cout << "Fbt_rho_mp11(): " << x << "   " << y << "   " << krho(x) << "   " << cos(k * x) << "   " << phimp11(x, y) << endl;
	}
	return ( 4.0 * i * exp(-i * k * y) * sum );
}

inline complex<double> Fbt_omega_pp00(double y, double k)
{
	double sum = 0.0;
	for (int ix = 0; ix < n_xi_pts; ix++)
	for (int iy = 0; iy < n_xi_pts; iy++)
	{
		double x = xi_pts_0_inf[ix];
		double y = xi_pts_0_inf[iy];
		sum += xi_wts_0_inf[ix] * xi_wts_0_inf[iy] * komega(x) * sin(k * x) * phipp00(x, y);
	}
	return ( 4.0 * i * exp(-i * k * y) * sum );
}

inline complex<double> Fbt_omega_pp20(double y, double k)
{
	double sum = 0.0;
	for (int ix = 0; ix < n_xi_pts; ix++)
	for (int iy = 0; iy < n_xi_pts; iy++)
	{
		double x = xi_pts_0_inf[ix];
		double y = xi_pts_0_inf[iy];
		sum += xi_wts_0_inf[ix] * xi_wts_0_inf[iy] * komega(x) * sin(k * x) * phipp20(x, y);
	}
	return ( 4.0 * i * exp(-i * k * y) * sum );
}

inline complex<double> Fbt_omega_pp02(double y, double k)
{
	double sum = 0.0;
	for (int ix = 0; ix < n_xi_pts; ix++)
	for (int iy = 0; iy < n_xi_pts; iy++)
	{
		double x = xi_pts_0_inf[ix];
		double y = xi_pts_0_inf[iy];
		sum += xi_wts_0_inf[ix] * xi_wts_0_inf[iy] * komega(x) * sin(k * x) * phipp02(x, y);
	}
	return ( 4.0 * i * exp(-i * k * y) * sum );
}

inline complex<double> Fbt_omega_mp11(double y, double k)
{
	double sum = 0.0;
	for (int ix = 0; ix < n_xi_pts; ix++)
	for (int iy = 0; iy < n_xi_pts; iy++)
	{
		double x = xi_pts_0_inf[ix];
		double y = xi_pts_0_inf[iy];
		sum += xi_wts_0_inf[ix] * xi_wts_0_inf[iy] * komega(x) * cos(k * x) * phimp11(x, y);
	}
	return ( 4.0 * exp(-i * k * y) * sum );
}

inline complex<double> Fbt_rho(double y, double k)
{
	double cy = cosh(y), sy = sinh(y);
	//phi0 = integrate_2D(phi0_int, xi_pts_minf_inf, xi_pts_minf_inf, xi_wts_minf_inf, xi_wts_minf_inf, n_xi_pts, n_xi_pts);
	phi0 = get_phi0(y);
	//cout << "Fbt_rho(): " << k << "   " << Fbt_rho_pp00(y, k) << "   " << Fbt_rho_pp20(y, k) << "   " << Fbt_rho_pp02(y, k) << "   " << Fbt_rho_mp11(y, k) << endl;
	return (
				Fbt_rho_pp00(y, k)
					+ cy*cy*Fbt_rho_pp20(y, k)
					+ sy*sy*Fbt_rho_pp02(y, k)
					+ sy*cy*Fbt_rho_mp11(y, k)
				//Fbt_rho_pp00(y, k)
				//	+ Fbt_rho_pp20(y, k)
				//	+ Fbt_rho_pp02(y, k)
				//	+ Fbt_rho_mp11(y, k)
			);
}

inline complex<double> Fbt_omega(double y, double k)
{
	double cy = cosh(y), sy = sinh(y);
	//phi0 = integrate_2D(phi0_int, xi_pts_minf_inf, xi_pts_minf_inf, xi_wts_minf_inf, xi_wts_minf_inf, n_xi_pts, n_xi_pts);
	phi0 = get_phi0(y);
	return (
				/*Fbt_omega_pp00(y, k)
					+ cy*cy*Fbt_omega_pp20(y, k)
					+ sy*sy*Fbt_omega_pp02(y, k)
					+ sy*cy*Fbt_omega_mp11(y, k)*/
				Fbt_omega_pp00(y, k)
					+ cy*cy*Fbt_omega_pp20(y, k)
					+ sy*sy*Fbt_omega_pp02(y, k)
					+ sy*cy*Fbt_omega_mp11(y, k)
			);
}

// End of file

#endif
