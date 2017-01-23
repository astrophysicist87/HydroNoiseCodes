#ifndef DEFS_H
#define DEFS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>
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

extern const double hbarC;

int integration_mode = 2;
bool include_baryon_chemical_potential_fluctations = true;
extern const double mu_pion;
extern double mu_part, mu_proton;

extern double vs, Neff, tauf, taui, Tf, Ti, nu, nuVB, ds, m, sf, s_at_mu_part;
extern double mByT, alpha0, psi0;
extern double chi_tilde_mu_mu, chi_tilde_T_mu, chi_tilde_T_T, Delta;

extern double screw_it_up_factor, current_kwt, current_DY;
extern int current_itau;

extern const int n_xi_pts;
extern const int n_k_pts;
extern const int n_tau_pts;
extern double * xi_pts_minf_inf, * xi_wts_minf_inf;
extern double * k_pts, * k_wts;
extern double * tau_pts, * tau_wts;
extern double * tau_pts_lower, * tau_wts_lower;
extern double * tau_pts_upper, * tau_wts_upper;
extern double * T_pts, * mu_pts;
extern double * T_pts_lower, * mu_pts_lower;
extern double * T_pts_upper, * mu_pts_upper;

extern double exp_delta, exp_gamma, exp_nu;
extern double T0, mu0, Tc, Pc, nc, sc, wc, muc;
extern double A0, A2, A4, C0, B, mui, muf, xi0, xibar0, etaBYs, RD, sPERn, Nf, qD, si, ni;
extern double a_at_tauf, vs2_at_tauf, vn2_at_tauf, vsigma2_at_tauf;

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

inline void set_phase_diagram_and_EOS_parameters()
{
	Nf = 2.0;											//number of massless flavors
	//T0 = 170.0;											//phase transition curve, T scale
	//mu0 = 1218.48;										//phase transition curve, mu scale
	T0 = 170.0 / hbarC;											//phase transition curve, T scale
	mu0 = 1218.48 / hbarC;										//phase transition curve, mu scale
	A4 = M_PI*M_PI*(16.0 + 10.5*Nf) / 90.0;				//coeff in P(T,mu) (i.e., EOS)
	A2 = Nf / 18.0;										//same
	A0 = Nf / (324.0 * M_PI * M_PI);					//same
	C0 = mu0*mu0*( A2 - 2.0*A0*mu0*mu0 / (T0*T0) );		//same
	B = 0.8 * pow(T0, 4.0);								//same

	//some parameters for the thermal conductivity critical enhancement
	//xi0 = 0.37;
	//qD = 1.0/xi0;
	//xibar0 = 0.69;
	qD = M_PI * T0;
	xi0 = 1.0 / qD;
	xibar0 = 0.701;
	RD = 1.05;
	exp_delta = 4.815;
	exp_gamma = 1.24;
	exp_nu = 0.63;
	etaBYs = 1.0 / (4.0 * M_PI);

	return;
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

inline void set_critical_point_parameters()
{
	//set critical point quantities
	//Tc = 160.0;
	//muc = 411.74;
	Tc = 160.0 / hbarC;
	muc = 411.74 / hbarC;
	Pc = P(Tc, muc);
	sc = s(Tc, muc);
	nc = n(Tc, muc);
	wc = w(Tc, muc);
	
	return;
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
		xibar0 * pow( t1+t2, -exp_nu/exp_gamma )		//fm^1
	);
}

inline double chi_B(double T, double mu)
{
	double t1 = (1.0/3.0)*((exp_delta - 1.0)/(2.0 - exp_gamma)) * pow(abs((T/Tc) - 1.0), exp_gamma);
	double t2 = 5.0 * exp_delta * pow(abs((n(T,mu)/nc) - 1.0), exp_delta - 1.0);
	return (
		( nc*nc / ( (exp_delta + 1.0) * Pc ) ) * pow(t1+t2, -1.0)		//MeV^2
	);
}

inline double cP(double T, double mu)
{
	double tp = (T/mu) * (2.0*A4*T*T + A2*mu*mu - C0) / (A2*T*T + 2.0*A0*mu*mu);
	double tm = 2.0*A2*T*mu / (A2*T*T + 6.0*A0*mu*mu);
	return (
		T * chi_B(T,mu) * (tp - tm) * (tp - tm)			//MeV^3
	);
}

inline double Delta_DT(double T, double mu)
{
	double xi_loc = xi(T, mu);
	return (
		RD * T * Omega(qD * xi_loc) / (6.0 * M_PI * etaBYs * s(T, mu) * xi_loc)
	);
}

inline double Delta_lambda(double T, double mu)
{
	return (cP(T,mu) * Delta_DT(T,mu));	//MeV^2
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

inline complex<double> beta(double k)
{
	return (
		std::sqrt( std::complex<double>( 0.25*(1.0-vsigma2_at_tauf)*(1.0-vsigma2_at_tauf) - vsigma2_at_tauf*k*k ) )
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
	
	double c1 = s_at_mu_part * chi_tilde_mu_mu;
	//double c2 = -sf * (chi_tilde_T_mu + chi_tilde_mu_mu * muf / Tf);
	double c2 = -s_at_mu_part * double(include_baryon_chemical_potential_fluctations) * (chi_tilde_T_mu + chi_tilde_mu_mu * mu_part / Tf);

	return ( (c1 * incompleteGamma4(mByT * cx) + c2 * incompleteGamma3(mByT * cx) ) / (cx*cx) );
}

inline double Fomega(double x)
{
	double cx = cosh(x);
	
	//return (Tf * tanh(x)*incompleteGamma4(mByT * cx) / (hbarC*cx*cx));
	return (Tf * tanh(x)*incompleteGamma4(mByT * cx) / (cx*cx));
}

inline double Fn(double x)
{
	double cx = cosh(x);
	
	double c1 = -s_at_mu_part * chi_tilde_T_mu;
	//double c2 = sf * (chi_tilde_T_T + chi_tilde_T_mu * muf / Tf);
	double c2 = s_at_mu_part * double(include_baryon_chemical_potential_fluctations) * (chi_tilde_T_T + chi_tilde_T_mu * mu_part / Tf);

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

inline complex<double> Gtilde_s(double k, double t_p)
{
	double T_loc = Tf;		//second argument always evaluated at t = tau_f
	double mu_loc = muf;	//second argument always evaluated at t = tau_f
	complex<double> b = beta(k);
	double t_by_tp = tauf / t_p;
	double f0 = vs2_at_tauf;
	double pref = (vsigma2_at_tauf-vs2_at_tauf)*pow(t_by_tp, -a_at_tauf);
	complex<double> f1 = (a_at_tauf/b)*sinh(b*log(t_by_tp)) + cosh(b*log(t_by_tp));

	return (
		(-i * k *mu_loc / (T_loc * vsigma2_at_tauf) ) * (f0 + pref*f1)	//dimensionless
	);
}

inline complex<double> Gtilde_omega(double k, double t_p)
{
	double T_loc = Tf;		//second argument always evaluated at t = tau_f
	double mu_loc = muf;	//second argument always evaluated at t = tau_f
	double t_by_tp = tauf / t_p;
	complex<double> b = beta(k);
	
	return (
		(k*k*sPERn/b) * ( vsigma2_at_tauf-vn2_at_tauf ) * pow(t_by_tp, -a_at_tauf) * sinh(b*log(t_by_tp))	//dimensionless
	);
}

inline complex<double> Gtilde_n(double k, double t_p)
{
	double T_loc = Tf;		//second argument always evaluated at t = tau_f
	double mu_loc = muf;	//second argument always evaluated at t = tau_f
	complex<double> b = beta(k);
	double t_by_tp = tauf / t_p;
	double f0 = vn2_at_tauf;
	double pref = (vsigma2_at_tauf-vn2_at_tauf)*pow(t_by_tp, -a_at_tauf);
	complex<double> f1 = (a_at_tauf/b)*sinh(b*log(t_by_tp)) + cosh(b*log(t_by_tp));

	/*if (current_itau == 45)
		cout << current_DY << "   " << k << "   " << current_kwt << "   " << t_p << "   " << tauf << "   " << T_loc << "   "
				<< mu_loc << "   " << a_at_tauf << "   " << b.real() << "   " << b.imag() << "   " << f0 << "  " << pref << "   "
				<< vsigma2_at_tauf << "   " << vn2_at_tauf << "   "
				<< ((a_at_tauf/b)*sinh(b*log(t_by_tp))).real() << "   " << ((a_at_tauf/b)*sinh(b*log(t_by_tp))).imag() << "   "
				<< cosh(b*log(t_by_tp)).real() << "   " << cosh(b*log(t_by_tp)).imag() << "   "
				<< k / vsigma2_at_tauf << endl;*/

	return (
		(i * k / vsigma2_at_tauf ) * (f0 + pref*f1)	 / (1.0 + screw_it_up_factor * k*k)		//dimensionless
	);
}

inline complex<double> tau_integration(complex<double> (*Gtilde_X)(double, double), complex<double> (*Gtilde_Y)(double, double), double k)
{
	complex<double> locsum(0,0);
	if (integration_mode == 1)
	{
		//lower section
		for (int it = 0; it < n_tau_pts; ++it)
		{
			double t_loc = tau_pts_lower[it];
			double T_loc = T_pts_lower[it];
			double mu_loc = mu_pts_lower[it];
			double arg = n(T_loc, mu_loc)*T_loc / ( s(T_loc, mu_loc)*w(T_loc, mu_loc) );
			locsum += tau_wts_lower[it] * pow(t_loc, -3.0) * Delta_lambda(T_loc, mu_loc) * pow(arg, 2.0)
					* (*Gtilde_X)(k, t_loc) * (*Gtilde_Y)(-k, t_loc);
		}
		//upper section
		for (int it = 0; it < n_tau_pts; ++it)
		{
			double t_loc = tau_pts_upper[it];
			double T_loc = T_pts_upper[it];
			double mu_loc = mu_pts_upper[it];
			double arg = n(T_loc, mu_loc)*T_loc / ( s(T_loc, mu_loc)*w(T_loc, mu_loc) );
			locsum += tau_wts_upper[it] * pow(t_loc, -3.0) * Delta_lambda(T_loc, mu_loc) * pow(arg, 2.0)
					* (*Gtilde_X)(k, t_loc) * (*Gtilde_Y)(-k, t_loc);
		}
	}
	else
	{
		//both sections together
		for (int it = 0; it < n_tau_pts; ++it)
		{
			double t_loc = tau_pts[it];
			double T_loc = T_pts[it];
			double mu_loc = mu_pts[it];
			double arg = n(T_loc, mu_loc)*T_loc / ( s(T_loc, mu_loc)*w(T_loc, mu_loc) );
			locsum += tau_wts[it] * pow(t_loc, -3.0) * Delta_lambda(T_loc, mu_loc) * pow(arg, 2.0)
					* (*Gtilde_X)(k, t_loc) * (*Gtilde_Y)(-k, t_loc);
		}
	}

	return (locsum);
}

inline complex<double> Ctilde_s_s(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_s, Gtilde_s, k);

	return ( sum );	//fm^2
}

inline complex<double> Ctilde_s_omega(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_s, Gtilde_omega, k);

	return ( sum );	//fm^2
}

inline complex<double> Ctilde_s_n(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_s, Gtilde_n, k);

	return ( sum );	//fm^2
}

inline complex<double> Ctilde_omega_s(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_omega, Gtilde_s, k);

	return ( sum );	//fm^2
}

inline complex<double> Ctilde_omega_omega(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_omega, Gtilde_omega, k);
	return ( sum );	//fm^2
}

inline complex<double> Ctilde_omega_n(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_omega, Gtilde_n, k);

	return ( sum );	//fm^2
}

inline complex<double> Ctilde_n_s(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_n, Gtilde_s, k);

	return ( sum );	//fm^2
}

inline complex<double> Ctilde_n_omega(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_n, Gtilde_omega, k);

	return ( sum );	//fm^2
}

inline complex<double> Ctilde_n_n(double k)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_n, Gtilde_n, k);

	return ( sum );	//fm^2
}

////////////////////////////////////////////////////////////////////////////////
// All functions in this section correspond to HBT correlator
////////////////////////////////////////////////////////////////////////////////

inline double zeta(double x1, double x2)
{
    return ( Tf / (cosh(x1) + cosh(x2)) );	//MeV^1
}

inline double gamma(double x1, double x2)
{
    double shx1_minus_shx2 = sinh(x1) - sinh(x2);
    return (shx1_minus_shx2 * shx1_minus_shx2);	//MeV^0
}

inline double f1(double x1, double x2)
{
	double z = zeta(x1, x2);				//MeV^1
	double gz = incompleteGamma4(m / z);	//MeV^0
	return (pow(z, 4.0) * gz);				//MeV^4
}

inline double f2(double x1, double x2)
{
	double z = zeta(x1, x2);																//MeV^1	
	double term1 = chi_tilde_mu_mu * cosh(x1) * pow(z, 5.0) * incompleteGamma5(m / z);		//MeV^3
	double term2 =  double(include_baryon_chemical_potential_fluctations)
					* (chi_tilde_mu_mu * mu_part + chi_tilde_T_mu * Tf)
					* pow(z, 4.0) * incompleteGamma4(m / z);								//MeV^3
	return (s_at_mu_part * (term1 - term2) / (Tf*Tf));												//MeV^4
}

inline double f3(double x1, double x2)
{
	double z = zeta(x1, x2);										//MeV^0
	return ( sinh(x1)*pow(z, 5.0)*incompleteGamma5(m / z) / Tf );	//MeV^4
}

inline double f4(double x1, double x2)
{
	double z = zeta(x1, x2);																//MeV^1
	double term1 = double(include_baryon_chemical_potential_fluctations)
					* (chi_tilde_T_mu*mu_part + chi_tilde_T_T*Tf)
						* pow(z, 4.0) * incompleteGamma4(m / z);							//MeV^3
	double term2 = chi_tilde_T_mu * cosh(x1) * pow(z, 5.0) * incompleteGamma5(m / z);		//MeV^3
	return (s_at_mu_part * (term1 - term2)/(Tf*Tf));													//MeV^4
}

inline double get_alpha0()
{
	double sum = 0.0;
	for (int ixi1 = 0; ixi1 < n_xi_pts; ++ixi1)
	for (int ixi2 = 0; ixi2 < n_xi_pts; ++ixi2)
	{
		double xi1 = xi_pts_minf_inf[ixi1];
		double xi2 = xi_pts_minf_inf[ixi2];
		sum += xi_wts_minf_inf[ixi1] * xi_wts_minf_inf[ixi2]
				* cosh(xi1) * cosh(xi2) * f1(xi1, xi2);			//MeV^4
	}
	return (1.0 / sum);											//MeV^{-4}
}

inline double get_psi0()
{
	double sum = 0.0;
	for (int ixi1 = 0; ixi1 < n_xi_pts; ++ixi1)
	for (int ixi2 = 0; ixi2 < n_xi_pts; ++ixi2)
	{
		double xi1 = xi_pts_minf_inf[ixi1];
		double xi2 = xi_pts_minf_inf[ixi2];
		sum += xi_wts_minf_inf[ixi1] * xi_wts_minf_inf[ixi2]
				* gamma(xi1, xi2) * cosh(xi1) * cosh(xi2) * f1(xi1, xi2);	//MeV^4
	}
	return (alpha0 * sum);													//MeV^0
}

inline complex<double> Fbtilde_s(double k)
{
	complex<double> sum(0,0);
	for (int ixi1 = 0; ixi1 < n_xi_pts; ++ixi1)
	for (int ixi2 = 0; ixi2 < n_xi_pts; ++ixi2)
	{
		double xi1 = xi_pts_minf_inf[ixi1];
		double xi2 = xi_pts_minf_inf[ixi2];
		sum += xi_wts_minf_inf[ixi1] * xi_wts_minf_inf[ixi2] * exp(-i * k * xi1)
				* cosh(xi1) * cosh(xi2) * f2(xi1, xi2) * (gamma(xi1, xi2) - psi0);	//MeV^4
	}

	return (alpha0*sum);															//MeV^0
}

inline complex<double> Fbtilde_omega(double k)
{
	complex<double> sum(0,0);
	for (int ixi1 = 0; ixi1 < n_xi_pts; ++ixi1)
	for (int ixi2 = 0; ixi2 < n_xi_pts; ++ixi2)
	{
		double xi1 = xi_pts_minf_inf[ixi1];
		double xi2 = xi_pts_minf_inf[ixi2];
		sum += xi_wts_minf_inf[ixi1] * xi_wts_minf_inf[ixi2] * exp(-i * k * xi1)
				* cosh(xi1) * cosh(xi2) * f3(xi1, xi2) * (gamma(xi1, xi2) - psi0);	//MeV^4
	}

	return (alpha0*sum);															//MeV^0
}

inline complex<double> Fbtilde_n(double k)
{
	complex<double> sum(0,0);
	for (int ixi1 = 0; ixi1 < n_xi_pts; ++ixi1)
	for (int ixi2 = 0; ixi2 < n_xi_pts; ++ixi2)
	{
		double xi1 = xi_pts_minf_inf[ixi1];
		double xi2 = xi_pts_minf_inf[ixi2];
		sum += xi_wts_minf_inf[ixi1] * xi_wts_minf_inf[ixi2] * exp(-i * k * xi1)
				* cosh(xi1) * cosh(xi2) * f4(xi1, xi2) * (gamma(xi1, xi2) - psi0);	//MeV^4
	}

	return (alpha0*sum);															//MeV^0
}

////////////////////////////////////////////////////////////////////////////////
// Functions/stuff to solve for time-dependence of T and mu and/or Tf and muf
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

//compute final time-step Tf and muf

int input_get_Tf_and_muf_f (const gsl_vector * x, void * params, gsl_vector * f)
{
	double tau_local = ((struct rparams *) params)->tau;

	const double x0 = gsl_vector_get (x, 0);	//T
	const double x1 = gsl_vector_get (x, 1);	//mu

	const double y0 = s(x0, x1)/n(x0, x1) - sPERn;	//defines fixed s/n curve
	const double y1 = P(x0,x1);						//defines P==0 curve

	gsl_vector_set (f, 0, y0);
	gsl_vector_set (f, 1, y1);

	return GSL_SUCCESS;
}

void compute_Tf_and_muf()
{
	const gsl_multiroot_fsolver_type *gsl_T;
	gsl_multiroot_fsolver *gsl_s;

	int status;
	size_t i, iter = 0;

	const size_t n = 2;
	struct rparams p = {taui};
	gsl_multiroot_function f = {&input_get_Tf_and_muf_f, n, &p};

	double x_init[2] = {Ti, mui};
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
	Tf = gsl_vector_get (gsl_s->x, 0);
	muf = gsl_vector_get (gsl_s->x, 1);

	gsl_multiroot_fsolver_free (gsl_s);
	gsl_vector_free (x);
	
	return;
}

//compute T and mu at each time step

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
	//both sections together
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
		//cout << "full: " << T_pts[it] << "   " << mu_pts[it] << endl;

		gsl_multiroot_fsolver_free (gsl_s);
		gsl_vector_free (x);
	}
	
	return;
}

void populate_T_and_mu_vs_tau_part2()
{
	//lower section
	for (int it = 0; it < n_tau_pts; ++it)
	{
		const gsl_multiroot_fsolver_type *gsl_T;
		gsl_multiroot_fsolver *gsl_s;

		int status;
		size_t i, iter = 0;

		const size_t n = 2;
		struct rparams p = {tau_pts_lower[it]};
		gsl_multiroot_function f = {&input_f, n, &p};

		double x_init[2] = {guess_T(tau_pts_lower[it]), guess_mu(tau_pts_lower[it])};
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
		T_pts_lower[it] = gsl_vector_get (gsl_s->x, 0);
		mu_pts_lower[it] = gsl_vector_get (gsl_s->x, 1);
		//cout << "lower: " << T_pts_lower[it] << "   " << mu_pts_lower[it] << endl;

		gsl_multiroot_fsolver_free (gsl_s);
		gsl_vector_free (x);
		//if (1) exit (0);
	}

	//upper section
	for (int it = 0; it < n_tau_pts; ++it)
	{
		const gsl_multiroot_fsolver_type *gsl_T;
		gsl_multiroot_fsolver *gsl_s;

		int status;
		size_t i, iter = 0;

		const size_t n = 2;
		struct rparams p = {tau_pts_upper[it]};
		gsl_multiroot_function f = {&input_f, n, &p};

		double x_init[2] = {guess_T(tau_pts_upper[it]), guess_mu(tau_pts_upper[it])};
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
		T_pts_upper[it] = gsl_vector_get (gsl_s->x, 0);
		mu_pts_upper[it] = gsl_vector_get (gsl_s->x, 1);
		//cout << "upper: " << T_pts_upper[it] << "   " << mu_pts_upper[it] << endl;

		gsl_multiroot_fsolver_free (gsl_s);
		gsl_vector_free (x);
	}
	
	return;
}

//a miscellaneous function to help split up tau integral
void break_up_integral(int nt, double &max_DL, double &tauc, double &width)
{
	double Delta_t = (tauf - taui - 1.0) / (double)nt;
	vector<double> Delta_lambda_vec;
	vector<double> tpts;

	//check Delta_lambda
	for (int it = 1; it < nt; ++it)
	{
		double t_loc = taui + 0.5 + (double)it * Delta_t;
		double T_loc = interpolate1D(tau_pts, T_pts, t_loc, n_tau_pts, 1, false);
		double mu_loc = interpolate1D(tau_pts, mu_pts, t_loc, n_tau_pts, 1, false);
		double DL = Delta_lambda(T_loc, mu_loc);
		Delta_lambda_vec.push_back(DL);
		tpts.push_back(t_loc);
		//cout << t_loc << "   " << DL << endl;
	}
	std::vector<double>::iterator result = std::max_element(Delta_lambda_vec.begin(), Delta_lambda_vec.end());
    int idx = std::distance(Delta_lambda_vec.begin(), result);
	max_DL = Delta_lambda_vec[idx];
	tauc = tpts[idx];
	width = 0.0;
	double den = 0.0;
	for (int ii = 0; ii < nt; ++ii)
	{
		width += pow(tpts[ii] - tauc, 2.0) * Delta_lambda_vec[ii];
		den += Delta_lambda_vec[ii];
	}
	width = sqrt(width / den);
	
	//cout << max_DL << "   " << tauc << "   " << width << endl;
	return;
}

// End of file

#endif
