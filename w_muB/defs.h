#ifndef DEFS_H
#define DEFS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>
#include <cmath>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

using namespace std;

#include "lib.h"

/*USAGE: debugger(__LINE__, __FILE__);*/
void inline debugger(int cln, const char* cfn)
{
	cerr << "You made it to " << cfn << ":" << cln << "!" << endl;
	return;
}

extern double vs, Neff, tauf, taui, Tf, Ti, nu, nuVB, ds, A, m, sf;
extern double mByT, alpha0, phi0;
extern double chi_tilde_mu_mu, chi_tilde_T_mu, chi_tilde_T_T, Delta;

extern const int n_xi_pts;
extern const int n_k_pts;
extern const int n_tau_pts;
extern double * xi_pts_0_inf, * xi_wts_0_inf, * xi_pts_minf_inf, * xi_wts_minf_inf;
extern double * k_pts, * k_wts;
extern double * tau_pts, * tau_wts;
extern double * T_pts, * mu_pts;

extern double exp_delta, exp_gamma, exp_nu;
extern double T0, mu0, Tc, Pc, nc, sc, wc, muc;
extern double A0, A2, A4, C0, B, mui, muf, xi0, xibar0, etaBYs, RD, sPERn, Nf, qD, si, ni;

//general functions
inline double Omega(double x)
{
	return (
		0.48*tanh(0.23*x) + (1.04 / M_PI) * atan(0.65*x)
	);
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

//equation of state and other thermodynamic relations
inline double P(double T, double mu)
{
	return(
		A4*pow(T, 4.0) + A2*T*T*mu*mu + A0*pow(mu, 4.0) - C0*T*T - B
	);
}

// functions to guess seed values for T and mu,
// followed by functions to iteratively solve equations
// for T and mu

inline double guess_T(double tau)
{
	return (Ti * pow(taui / tau, 1.0/3.0));
}

inline double guess_mu(double tau)
{
	return (mui * pow(taui / tau, 1.0/3.0));
}

///////////////////////////////////////////////////////////
//two separate definitions of s and n each, for convenience
///////////////////////////////////////////////////////////
inline double s(double tau)
{
	return (si * taui / tau);
}

inline double n(double tau)
{
	return (ni * taui / tau);
}

inline double s(double T, double mu)
{
	return (
		-2.0 * C0 * T + 4.0 * A4 * pow(T, 3.0) + 2.0 * A2 * T * mu * mu
	);
}

inline double n(double T, double mu)
{
	return (
		2.0 * A2 * T * T * mu + 4.0 * A0 * pow(mu, 3.0)
	);
}
///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

inline double w(double T, double mu)
{
	return (T * s(T, mu) + mu * n(T, mu));
}

inline double chi_TT(double T, double mu)
{
	return (
		2.0 * ( -C0 + 6.0 * A4 * T * T + A2 * mu * mu )
	);
}

inline double chi_Tmu(double T, double mu)
{
	return (
		4.0 * A2 * T * mu
	);
}

inline double chi_mumu(double T, double mu)
{
	return (
		2.0 * ( A2 * T * T + 6.0 * A0 * mu * mu )
	);
}


inline double xi(double T, double mu)
{
	double t1 = (1.0/3.0)*((exp_delta - 1.0)/(2.0 - exp_gamma)) * pow(abs((T/Tc) - 1.0), exp_gamma);
	double t2 = 5.0 * exp_delta * pow(abs((n(T,mu)/nc) - 1.0), exp_delta - 1.0);
	return (
		xibar0 * pow( t1+t2, -exp_nu/exp_gamma )
	);
}

inline double chi_B(double T, double mu)
{
	double t1 = (1.0/3.0)*((exp_delta - 1.0)/(2.0 - exp_gamma)) * pow(abs((T/Tc) - 1.0), exp_gamma);
	double t2 = 5.0 * exp_delta * pow(abs((n(T,mu)/nc) - 1.0), exp_delta - 1.0);
	return (
		(nc*nc / ((exp_delta + 1.0) * Pc)) * pow(t1+t2, -1.0)
	);
}

inline double cV(double T, double mu)
{
	double cTT = chi_TT(T, mu);
	double cTmu = chi_Tmu(T, mu);
	double cmumu = chi_mumu(T, mu);

	return (
		T * (cTT*cmumu - cTmu*cTmu) / cmumu
	);
}

inline double cP(double T, double mu)
{
	double tp = (T/mu) * (2.0*A4*T*T + A2*mu*mu - C0) / (A2*T*T + 2.0*A0*mu*mu);
	double tm = 2.0*A2*T*mu / (A2*T*T + 6.0*A0*mu*mu);
	return (
		0.0*cV(T,mu) + T * chi_B(T,mu) * (tp - tm) * (tp - tm)
	);
}
inline double Delta_DT(double T, double mu)
{
	double xi_loc = xi(T, mu);
	return (
		197.33 * RD * T * Omega(qD * xi_loc) / (6.0 * M_PI * etaBYs * s(T, mu) * xi_loc)
	);
}

inline double Delta_lambda(double T, double mu)
{
	return (cP(T,mu) * Delta_DT(T,mu));
}

//speeds of sound for FTd-Green functions
inline double d(double T, double mu)
{
	return (
		2.0*A4*A2*pow(T, 4.0) + (12.0*A4*A0 - A2*A2)*mu*mu*T*T + 2.0*A2*A0*pow(mu, 4.0)
			- C0*(A2*T*T + 6.0*A0*mu*mu) / 3.0
	);
}

inline double vn2(double T, double mu)
{
	return (
		(1.0/3.0) - ( 2.0*C0*( A2*T*T + 6.0*A0*mu*mu ) ) / ( 9.0*d(T, mu) )
	);
}

inline double vs2(double T, double mu)
{
	return (
		(1.0/3.0) + ( ( 4.0*A2*C0*T*T ) / ( 9.0*d(T, mu) ) )
	);
}

inline double vsigma2(double T, double mu)
{
	double s_loc = s(T, mu);
	double n_loc = n(T, mu);
	double w_loc = w(T, mu);
	double vn2_loc = vn2(T, mu);
	double vs2_loc = vs2(T, mu);
	return (
		(T*s_loc*vn2_loc + mu*n_loc*vs2_loc) / w_loc
	);
}















inline double alpha(double T, double mu)
{
	return (
		0.5 * (1.0 - vsigma2(T, mu))
	);
}

inline complex<double> beta(double T, double mu, double k)
{
	double vsig2 = vsigma2(T, mu);
	return (
		std::sqrt( std::complex<double>( 0.25*(1.0-vsig2)*(1.0-vsig2) - vsig2*k*k ) )
	);
}

inline double norm_int(double x)
{
	double cx = cosh(x);
	return (incompleteGamma3(mByT * cx) / (cx*cx));
}

inline double Fs(double x)
{
	double cx = cosh(x);
	
	double c1 = sf * chi_tilde_mu_mu;
	double c2 = -sf * (chi_tilde_T_mu + chi_tilde_mu_mu * muf / Tf);

	return ( (c1 * incompleteGamma4(mByT * cx) + c2 * incompleteGamma3(mByT * cx) ) / (cx*cx) );
}

inline double Fomega(double x)
{
	double cx = cosh(x);
	
	return (Tf * tanh(x)*incompleteGamma4(mByT * cx) / (cx*cx));
}

inline double Fn(double x)
{
	double cx = cosh(x);
	
	double c1 = -sf * chi_tilde_T_mu;
	double c2 = sf * (chi_tilde_T_T + chi_tilde_T_mu * muf / Tf);

	return ( (c1 * incompleteGamma4(mByT * cx) + c2 * incompleteGamma3(mByT * cx) ) / (cx*cx) );
}

inline complex<double> Ftilde_s(double k)
{
	return (integrate_1D_FT(Fs, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts, k));
}

inline complex<double> Ftilde_omega(double k)
{
	return (integrate_1D_FT(Fomega, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts, k));
}

inline complex<double> Ftilde_n(double k)
{
	return (integrate_1D_FT(Fn, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts, k));
}

//inline complex<double> Gtilde_s(double k, double t, double t_p)
inline complex<double> Gtilde_s(double k, double t, int i_t_p)
{
	//double T_loc = T(t);
	//double mu_loc = mu(t);
	double t_p = tau_pts[i_t_p];
	double T_loc = Tf;		//second argument always evaluated at t = tau_f
	double mu_loc = muf;	//second argument always evaluated at t = tau_f
	double a = alpha(T_loc, mu_loc);
	complex<double> b = beta(T_loc, mu_loc, k);
	double t_by_tp = t / t_p;
	double vs2_at_t = vs2(T_loc, mu_loc);
	double vsigma2_at_t = vsigma2(T_loc, mu_loc);
	double f0 = vs2_at_t;
	double pref = (vsigma2_at_t-vs2_at_t)*pow(t_by_tp, -a);
	complex<double> f1 = (a/b)*sinh(b*log(t_by_tp)) + cosh(b*log(t_by_tp));

	return (
		(-i * k *mu_loc / (T_loc * vsigma2_at_t) ) * (f0 + pref*f1)
	);
}

//inline complex<double> Gtilde_omega(double k, double t, double t_p)
inline complex<double> Gtilde_omega(double k, double t, int i_t_p)
{
	//double T_loc = T(t);
	//double mu_loc = mu(t);
	double t_p = tau_pts[i_t_p];
	double T_loc = Tf;		//second argument always evaluated at t = tau_f
	double mu_loc = muf;	//second argument always evaluated at t = tau_f
	double a = alpha(T_loc, mu_loc);
	double t_by_tp = t / t_p;
	complex<double> b = beta(T_loc, mu_loc, k);
	double vn2_at_t = vn2(T_loc, mu_loc);
	double vsigma2_at_t = vsigma2(T_loc, mu_loc);
	
	return (
		(k*k*sPERn/b) * ( vsigma2_at_t-vn2_at_t ) * pow(t_by_tp, -a) * sinh(b*log(t_by_tp))
	);
}

//inline complex<double> Gtilde_n(double k, double t, double t_p)
inline complex<double> Gtilde_n(double k, double t, int i_t_p)
{
	//double T_loc = T(t);
	//double mu_loc = mu(t);
	double t_p = tau_pts[i_t_p];
	double T_loc = Tf;		//second argument always evaluated at t = tau_f
	double mu_loc = muf;	//second argument always evaluated at t = tau_f
	double a = alpha(T_loc, mu_loc);
	complex<double> b = beta(T_loc, mu_loc, k);
	double t_by_tp = t / t_p;
	double vn2_at_t = vn2(T_loc, mu_loc);
	double vsigma2_at_t = vsigma2(T_loc, mu_loc);
	double f0 = vn2_at_t;
	double pref = (vsigma2_at_t-vn2_at_t)*pow(t_by_tp, -a);
	complex<double> f1 = (a/b)*sinh(b*log(t_by_tp)) + cosh(b*log(t_by_tp));

	return (
		(i * k / vsigma2_at_t ) * (f0 + pref*f1)
	);
}

inline complex<double> Ctilde_s_s(double k)
{
	complex<double> sum(0,0);
	for (int it = 0; it < n_tau_pts; ++it)
	{
		double t_loc = tau_pts[it];
		//double T_loc = T(t_loc);
		//double mu_loc = mu(t_loc);
		double T_loc = T_pts[it];
		double mu_loc = mu_pts[it];
		if (abs( 1.0 - ni*taui/t_loc/n(T_loc, mu_loc) ) > 1.e-6 || abs(1.0 - si*taui/t_loc/s(T_loc, mu_loc) ) > 1.e-6 )
			cout << n(T_loc, mu_loc) << "   " << ni*taui/t_loc << "   " << s(T_loc, mu_loc) << "   " << si*taui/t_loc << endl;
		sum += tau_wts[it] * pow(t_loc, -3.0) * Delta_lambda(T_loc, mu_loc) * pow(n(T_loc, mu_loc)*T_loc / ( s(T_loc, mu_loc)*w(T_loc, mu_loc) ), 2.0)
				* Gtilde_s(k, tauf, it) * Gtilde_s(-k, tauf, it);
				//* Gtilde_s(k, tauf, t_loc) * Gtilde_s(-k, tauf, t_loc);
	}

	return (
		sum/A
	);
}

inline complex<double> Ctilde_s_omega(double k)
{
	complex<double> sum(0,0);
	for (int it = 0; it < n_tau_pts; ++it)
	{
		double t_loc = tau_pts[it];
		//double T_loc = T(t_loc);
		//double mu_loc = mu(t_loc);
		double T_loc = T_pts[it];
		double mu_loc = mu_pts[it];
		if (abs( n(T_loc, mu_loc) - ni*taui/t_loc ) > 1.e-6 || abs( s(T_loc, mu_loc) - si*taui/t_loc ) > 1.e-6 )
			cout << n(T_loc, mu_loc) << "   " << ni*taui/t_loc << "   " << s(T_loc, mu_loc) << "   " << si*taui/t_loc << endl;
		sum += tau_wts[it] * pow(t_loc, -3.0) * Delta_lambda(T_loc, mu_loc) * pow(n(T_loc, mu_loc)*T_loc / ( s(T_loc, mu_loc)*w(T_loc, mu_loc) ), 2.0)
				* Gtilde_s(k, tauf, it) * Gtilde_omega(-k, tauf, it);
				//* Gtilde_s(k, tauf, t_loc) * Gtilde_omega(-k, tauf, t_loc);
	}

	return (
		sum/A
	);
}

inline complex<double> Ctilde_s_n(double k)
{
	complex<double> sum(0,0);
	for (int it = 0; it < n_tau_pts; ++it)
	{
		double t_loc = tau_pts[it];
		//double T_loc = T(t_loc);
		//double mu_loc = mu(t_loc);
		double T_loc = T_pts[it];
		double mu_loc = mu_pts[it];
		if (abs( n(T_loc, mu_loc) - ni*taui/t_loc ) > 1.e-6 || abs( s(T_loc, mu_loc) - si*taui/t_loc ) > 1.e-6 )
			cout << n(T_loc, mu_loc) << "   " << ni*taui/t_loc << "   " << s(T_loc, mu_loc) << "   " << si*taui/t_loc << endl;
		sum += tau_wts[it] * pow(t_loc, -3.0) * Delta_lambda(T_loc, mu_loc) * pow(n(T_loc, mu_loc)*T_loc / ( s(T_loc, mu_loc)*w(T_loc, mu_loc) ), 2.0)
				* Gtilde_s(k, tauf, it) * Gtilde_n(-k, tauf, it);
				//* Gtilde_s(k, tauf, t_loc) * Gtilde_n(-k, tauf, t_loc);
	}

	return (
		sum/A
	);
}

inline complex<double> Ctilde_omega_s(double k)
{
	complex<double> sum(0,0);
	for (int it = 0; it < n_tau_pts; ++it)
	{
		double t_loc = tau_pts[it];
		//double T_loc = T(t_loc);
		//double mu_loc = mu(t_loc);
		double T_loc = T_pts[it];
		double mu_loc = mu_pts[it];
		if (abs( n(T_loc, mu_loc) - ni*taui/t_loc ) > 1.e-6 || abs( s(T_loc, mu_loc) - si*taui/t_loc ) > 1.e-6 )
			cout << n(T_loc, mu_loc) << "   " << ni*taui/t_loc << "   " << s(T_loc, mu_loc) << "   " << si*taui/t_loc << endl;
		sum += tau_wts[it] * pow(t_loc, -3.0) * Delta_lambda(T_loc, mu_loc) * pow(n(T_loc, mu_loc)*T_loc / ( s(T_loc, mu_loc)*w(T_loc, mu_loc) ), 2.0)
				* Gtilde_omega(k, tauf, it) * Gtilde_s(-k, tauf, it);
				//* Gtilde_omega(k, tauf, t_loc) * Gtilde_s(-k, tauf, t_loc);
	}

	return (
		sum/A
	);
}

inline complex<double> Ctilde_omega_omega(double k)
{
	complex<double> sum(0,0);
	for (int it = 0; it < n_tau_pts; ++it)
	{
		double t_loc = tau_pts[it];
		//double T_loc = T(t_loc);
		//double mu_loc = mu(t_loc);
		double T_loc = T_pts[it];
		double mu_loc = mu_pts[it];
		if (abs( n(T_loc, mu_loc) - ni*taui/t_loc ) > 1.e-6 || abs( s(T_loc, mu_loc) - si*taui/t_loc ) > 1.e-6 )
			cout << n(T_loc, mu_loc) << "   " << ni*taui/t_loc << "   " << s(T_loc, mu_loc) << "   " << si*taui/t_loc << endl;
		sum += tau_wts[it] * pow(t_loc, -3.0) * Delta_lambda(T_loc, mu_loc) * pow(n(T_loc, mu_loc)*T_loc / ( s(T_loc, mu_loc)*w(T_loc, mu_loc) ), 2.0)
				* Gtilde_omega(k, tauf, it) * Gtilde_omega(-k, tauf, it);
				//* Gtilde_omega(k, tauf, t_loc) * Gtilde_omega(-k, tauf, t_loc);
	}

	return (
		sum/A
	);
}

inline complex<double> Ctilde_omega_n(double k)
{
	complex<double> sum(0,0);
	for (int it = 0; it < n_tau_pts; ++it)
	{
		double t_loc = tau_pts[it];
		//double T_loc = T(t_loc);
		//double mu_loc = mu(t_loc);
		double T_loc = T_pts[it];
		double mu_loc = mu_pts[it];
		if (abs( n(T_loc, mu_loc) - ni*taui/t_loc ) > 1.e-6 || abs( s(T_loc, mu_loc) - si*taui/t_loc ) > 1.e-6 )
			cout << n(T_loc, mu_loc) << "   " << ni*taui/t_loc << "   " << s(T_loc, mu_loc) << "   " << si*taui/t_loc << endl;
		sum += tau_wts[it] * pow(t_loc, -3.0) * Delta_lambda(T_loc, mu_loc) * pow(n(T_loc, mu_loc)*T_loc / ( s(T_loc, mu_loc)*w(T_loc, mu_loc) ), 2.0)
				* Gtilde_omega(k, tauf, it) * Gtilde_n(-k, tauf, it);
				//* Gtilde_omega(k, tauf, t_loc) * Gtilde_n(-k, tauf, t_loc);
	}

	return (
		sum/A
	);
}

inline complex<double> Ctilde_n_s(double k)
{
	complex<double> sum(0,0);
	for (int it = 0; it < n_tau_pts; ++it)
	{
		double t_loc = tau_pts[it];
		//double T_loc = T(t_loc);
		//double mu_loc = mu(t_loc);
		double T_loc = T_pts[it];
		double mu_loc = mu_pts[it];
		if (abs( n(T_loc, mu_loc) - ni*taui/t_loc ) > 1.e-6 || abs( s(T_loc, mu_loc) - si*taui/t_loc ) > 1.e-6 )
			cout << n(T_loc, mu_loc) << "   " << ni*taui/t_loc << "   " << s(T_loc, mu_loc) << "   " << si*taui/t_loc << endl;
		sum += tau_wts[it] * pow(t_loc, -3.0) * Delta_lambda(T_loc, mu_loc) * pow(n(T_loc, mu_loc)*T_loc / ( s(T_loc, mu_loc)*w(T_loc, mu_loc) ), 2.0)
				* Gtilde_n(k, tauf, it) * Gtilde_s(-k, tauf, it);
				//* Gtilde_n(k, tauf, t_loc) * Gtilde_s(-k, tauf, t_loc);
	}

	return (
		sum/A
	);
}

inline complex<double> Ctilde_n_omega(double k)
{
	complex<double> sum(0,0);
	for (int it = 0; it < n_tau_pts; ++it)
	{
		double t_loc = tau_pts[it];
		//double T_loc = T(t_loc);
		//double mu_loc = mu(t_loc);
		double T_loc = T_pts[it];
		double mu_loc = mu_pts[it];
		if (abs( n(T_loc, mu_loc) - ni*taui/t_loc ) > 1.e-6 || abs( s(T_loc, mu_loc) - si*taui/t_loc ) > 1.e-6 )
			cout << n(T_loc, mu_loc) << "   " << ni*taui/t_loc << "   " << s(T_loc, mu_loc) << "   " << si*taui/t_loc << endl;
		sum += tau_wts[it] * pow(t_loc, -3.0) * Delta_lambda(T_loc, mu_loc) * pow(n(T_loc, mu_loc)*T_loc / ( s(T_loc, mu_loc)*w(T_loc, mu_loc) ), 2.0)
				* Gtilde_n(k, tauf, it) * Gtilde_omega(-k, tauf, it);
				//* Gtilde_n(k, tauf, t_loc) * Gtilde_omega(-k, tauf, t_loc);
	}

	return (
		sum/A
	);
}

inline complex<double> Ctilde_n_n(double k)
{
	complex<double> sum(0,0);
	for (int it = 0; it < n_tau_pts; ++it)
	{
		double t_loc = tau_pts[it];
		//double T_loc = T(t_loc);
		//double mu_loc = mu(t_loc);
		double T_loc = T_pts[it];
		double mu_loc = mu_pts[it];
		if (abs( n(T_loc, mu_loc) - ni*taui/t_loc ) > 1.e-6 || abs( s(T_loc, mu_loc) - si*taui/t_loc ) > 1.e-6 )
			cout << n(T_loc, mu_loc) << "   " << ni*taui/t_loc << "   " << s(T_loc, mu_loc) << "   " << si*taui/t_loc << endl;
		sum += tau_wts[it] * pow(t_loc, -3.0) * Delta_lambda(T_loc, mu_loc) * pow(n(T_loc, mu_loc)*T_loc / ( s(T_loc, mu_loc)*w(T_loc, mu_loc) ), 2.0)
				* Gtilde_n(k, tauf, it) * Gtilde_n(-k, tauf, it);
				//* Gtilde_n(k, tauf, t_loc) * Gtilde_n(-k, tauf, t_loc);
	}

	return (
		sum/A
	);
}

////////////////////////////////////////////////////////////////////////////////
// All functions in this section correspond to HBT correlator
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


////////////////////////////////////////////////////////////////////////////////
// Functions/stuff to solve for time-dependence of T and mu
////////////////////////////////////////////////////////////////////////////////

struct rparams
{
	double tau;
};

int print_state (size_t iter, gsl_multiroot_fsolver * s)
{
	printf ("iter = %3u x = % .3f % .3f "
		"f(x) = % .3e % .3e\n",
		iter,
		gsl_vector_get (s->x, 0), 
		gsl_vector_get (s->x, 1),
		gsl_vector_get (s->f, 0), 
		gsl_vector_get (s->f, 1));
}

int input_f (const gsl_vector * x, void * params, gsl_vector * f)
{
	double tau_local = ((struct rparams *) params)->tau;

	const double x0 = gsl_vector_get (x, 0);	//T
	const double x1 = gsl_vector_get (x, 1);	//mu

	const double y0 = s(x0, x1) - s(tau_local);
	const double y1 = n(x0, x1) - n(tau_local);

	gsl_vector_set (f, 0, y0);
	gsl_vector_set (f, 1, y1);

	return GSL_SUCCESS;
}

void populate_T_and_mu_vs_tau()
{
	for (int it = 0; it < n_tau_pts; ++it)
	{
		const gsl_multiroot_fsolver_type *gsl_T;
		gsl_multiroot_fsolver *gsl_s;

		int status;
		size_t i, iter = 0;

		const size_t n = 2;
		struct rparams p = {tau_pts[it]};
		gsl_multiroot_function f = {&input_f, n, &p};

		double x_init[2] = {guess_T(tau_pts[it]), guess_mu(tau_pts[it])};
		gsl_vector *x = gsl_vector_alloc (n);

		gsl_vector_set (x, 0, x_init[0]);
		gsl_vector_set (x, 1, x_init[1]);

		gsl_T = gsl_multiroot_fsolver_hybrids;
		gsl_s = gsl_multiroot_fsolver_alloc (gsl_T, 2);
		gsl_multiroot_fsolver_set (gsl_s, &f, x);

		//print_state (iter, gsl_s);

		do
		{
			iter++;
			status = gsl_multiroot_fsolver_iterate (gsl_s);

			//print_state (iter, gsl_s);

			if (status)   /* check if solver is stuck */
				break;

			status = gsl_multiroot_test_residual (gsl_s->f, 1e-7);
		}
		while (status == GSL_CONTINUE && iter < 1000);

		//printf ("status = %s\n", gsl_strerror (status));

		//finally, store results
		T_pts[it] = gsl_vector_get (gsl_s->x, 0);
		mu_pts[it] = gsl_vector_get (gsl_s->x, 1);

		gsl_multiroot_fsolver_free (gsl_s);
		gsl_vector_free (x);
	}
	
	return;
}

// End of file

#endif
