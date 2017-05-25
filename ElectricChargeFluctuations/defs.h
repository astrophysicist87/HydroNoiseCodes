#ifndef DEFS_H
#define DEFS_H

#include <iostream>
#include <fstream>
#include <sstream>
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

string truestring = "true";
string falsestring = "false";

inline string return_boolean_string(bool test){return (test ? truestring : falsestring);}

struct chosen_particle
{
	int index;
	double mass;
};

extern const double hbarC;
extern const double Cem;

double fraction_of_evolution;

//bool use_lattice_transport = false;

//double two_pi_DQ_T = 1.0;	//default constant value

extern long n_interp;

extern double vs, Neff, tauf, taui, Tf, Ti, nu, nuVB, ds, sf;
extern double chi_tilde_mu_mu, chi_tilde_T_mu, chi_tilde_T_T, Delta;

extern double chi_T_T, chi_T_mu, chi_mu_mu;

extern const int n_xi_pts;
extern const int n_k_pts;
extern const int n_tau_pts;
extern double * xi_pts_minf_inf, * xi_wts_minf_inf;
extern double * k_pts, * k_wts;
extern double * tau_pts, * tau_wts;
extern double * T_pts;

chosen_particle particle1;
chosen_particle particle2;

int current_ik = -1;

extern double * interp_T_pts, * interp_transport_pts;

extern double T0, mu0, Tc, Pc, nc, sc, wc, muc;
extern double A0, A2, A4, C0, B, xi0, xibar0, etaBYs, RD, sPERn, Nf, qD, si;

//general functions
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

long get_filelength(string filename)
{
   long length=0; 
   char line[512];
   ostringstream filepath_stream;
   filepath_stream << "./" << filename;
   ifstream infile(filepath_stream.str().c_str());
   //Determine the length of the file
   while (!infile.eof ())
   {
      infile.getline(line, 512);
      length++;
   }
   length = length-1;
   infile.close();
   return(length);
}

inline void read_in_transport_table(long filelength, string input_filename)
{
	ostringstream filename_stream;
	filename_stream << "./" << input_filename;
	ifstream in(filename_stream.str().c_str());

	for (int iline = 0; iline < filelength; ++iline)
	{
		in >> interp_T_pts[iline] >> interp_transport_pts[iline];
		interp_T_pts[iline] /= hbarC;
	}

	in.close();
	return;
}

//equation of state and other thermodynamic relations
inline double P(double T)
{
	return(
		A4*pow(T, 4.0) - C0*T*T - B
	);
}

inline void set_phase_diagram_and_EOS_parameters()
{
	Nf = 3.0;											//number of massless flavors
	T0 = 170.0 / hbarC;											//phase transition curve, T scale
	mu0 = 1218.48 / hbarC;										//phase transition curve, mu scale
	A4 = M_PI*M_PI*(16.0 + 10.5*Nf) / 90.0;				//coeff in P(T) (i.e., EOS)
	A2 = Nf / 18.0;										//same
	A0 = Nf / (324.0 * M_PI * M_PI);					//same
	C0 = mu0*mu0*( A2 - 2.0*A0*mu0*mu0 / (T0*T0) );		//same
	B = 0.8 * pow(T0, 4.0);								//same

	return;
}

// functions to guess seed value for T,
// followed by functions to iteratively solve equations
// for T

inline double guess_T(double tau)
{
	return (Ti * pow(taui / tau, 1.0/3.0));
}

///////////////////////////////////////////////////////////
//two separate definitions of s, for convenience
///////////////////////////////////////////////////////////
inline double s_vs_tau(double tau)
{
	return (si * taui / tau);
}

inline double s_vs_T(double T)
{
	return (
		-2.0 * C0 * T + 4.0 * A4 * T*T*T
	);
}

///////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////

inline double w(double T)
{
	return (T * s_vs_T(T));
}

inline double chi_TT(double T)
{
	return ( 2.0 * ( -C0 + 6.0 * A4 * T * T ) );
}

inline double chi_Tmu(double T)
{
	T;
	return ( 0.0 );
}

inline double chi_mumu(double T)
{
	return ( 4.0 * A2 * T * T );
}

inline double norm_int(double x, void * p)
{
	struct chosen_particle * params = (struct chosen_particle *)p;
	double cx = cosh(x);
	double mByT = (params->mass) / Tf;
	return (incompleteGamma3(mByT * cx) / (cx*cx));
}

inline double Fn(double x, void * p)
{
	struct chosen_particle * params = (struct chosen_particle *)p;
	double cx = cosh(x);
	
	double c1 = 0.0;
	double c2 = sf * chi_tilde_T_T;

	double mByT = (params->mass) / Tf;

	return ( (c1 * incompleteGamma4(mByT * cx) + c2 * incompleteGamma3(mByT * cx) ) / (cx*cx) );
}

inline complex<double> Ftilde_n(double k, void * p)
{
	return (integrate_1D_FT(Fn, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts, k, p));
}

inline complex<double> Gtilde_n(double k)
{
	return ( i * k ); //dimensionless
}

inline complex<double> tau_integration(complex<double> (*Gtilde_X)(double), complex<double> (*Gtilde_Y)(double), double k, bool use_lattice_transport, double two_pi_DQ_T)
{
	complex<double> result(0,0);
	if (use_lattice_transport)
	{
		for (int it = 0; it < n_tau_pts; ++it)
		{
			double tau = tau_pts[it];
			double T_loc = T_pts[it];
			double s_loc = s_vs_T(T_loc);
			double tau3 = tau*tau*tau;
			double sigma_By_Cem_T = interpolate1D(interp_T_pts, interp_transport_pts, T_loc, n_interp, 0, 1);
			result += 4.0 * M_PI * Cem * tau_wts[it] * T_loc*T_loc*sigma_By_Cem_T / (s_loc*s_loc*tau3)
						* (*Gtilde_X)(k) * (*Gtilde_Y)(-k);
		}
	}
	else	//just use a constant value of 2*pi*D_Q*T and chi_Q
	{
		for (int it = 0; it < n_tau_pts; ++it)
		{
			double tau = tau_pts[it];
			double T_loc = T_pts[it];
			double s_loc = s_vs_T(T_loc);
			double tau3 = tau*tau*tau;
			double chi_Q = chi_mumu(T_loc);
			complex<double> tmp_result = 2.0 * tau_wts[it] * two_pi_DQ_T * chi_Q / (s_loc*s_loc*tau3)
						* (*Gtilde_X)(k) * (*Gtilde_Y)(-k);
//if (current_ik==0 && abs(two_pi_DQ_T - 1.0) < 1.e-10) cerr << "CHECK: " << tau << "   " << T_loc << "   " << s_loc << "   " << chi_Q << "   " << tmp_result.real() << endl;
			result += tmp_result;
		}
	}

	return (result);
}

inline void set_running_transport_integral(double * running_integral_array)
{
	const int n_x_pts = 11;
	double * x_pts = new double [n_x_pts];
	double * x_wts = new double [n_x_pts];
	gauss_quadrature(n_x_pts, 1, 0.0, 0.0, -1.0, 1.0, x_pts, x_wts);

	running_integral_array[0] = 0.0;
	for (int it = 0; it < 2*n_tau_pts-1; ++it)
	{
		double sum = 0.0;
		double t0 = all_tau_pts[it], t1 = all_tau_pts[it+1];
		double cen = 0.5*(t0+t1);
		double hw = 0.5*(t1-t0);
		for (int ix = 0; ix < n_x_pts; ++ix)
		{
			double t_loc = cen + hw * x_pts[ix];
			double T_loc = interpolate1D(all_tau_pts, all_T_pts, t_loc, 2*n_tau_pts, 0, false);
			double mu_loc = interpolate1D(all_tau_pts, all_mu_pts, t_loc, 2*n_tau_pts, 0, false);
			double arg = n(T_loc, mu_loc)*T_loc / ( s(T_loc, mu_loc)*w(T_loc, mu_loc) );
			sum += x_wts[ix] * hw * exp(t_loc / tauC) * Delta_lambda(T_loc, mu_loc) * pow(arg, 2.0) / t_loc;
			//sum += x_wts[ix] * hw * exp(t_loc / tauC);	//use this line to check that integration works correctly
		}
		running_integral_array[it+1] = exp(-t1 / tauC) * ( tauC * exp(t0 / tauC) * running_integral_array[it] + sum ) / tauC;	//this array contains eta(x) (defined on my whiteboard)
	}

	delete [] x_pts;
	delete [] x_wts;

	return;
}

inline complex<double> colored_tau_integration(complex<double> (*Gtilde_X)(double, double), complex<double> (*Gtilde_Y)(double, double), double k)
{
	complex<double> locsum(0,0);
	
	const int n_x_pts = 201;	//try this
	double * x_pts = new double [n_x_pts];
	double * x_wts = new double [n_x_pts];
	gauss_quadrature(n_x_pts, 1, 0.0, 0.0, -1.0, 1.0, x_pts, x_wts);

	double delta_tau_lower = -10.0 * tauC, delta_tau_upper = 10.0 * tauC;	//this bounds the interval where the integrand is large-ish
	for (int itp = 0; itp < 2*n_tau_pts; ++itp)
	{
		double tX_loc = all_tau_pts[itp];
		double TX_loc = all_T_pts[itp];
		double muX_loc = all_mu_pts[itp];
		complex<double> factor_X = (*Gtilde_X)(k, tX_loc) / tX_loc;

		double tau_lower = max(taui, tX_loc + delta_tau_lower);			//if lower limit goes before beginning of lifetime, just start at tau0
		double tau_upper = min(tauf, tX_loc + delta_tau_upper);			//if upper limit goes past end of lifetime, just end at tauf
		double hw_loc = 0.5 * (tau_upper - tau_lower);
		double cen_loc = 0.5 * (tau_upper + tau_lower);

		complex<double> sum_X = 0.0;
		for (int ix = 0; ix < n_x_pts; ++ix)
		{
			double tY_loc = cen_loc + hw_loc * x_pts[ix];
			double TY_loc = interpolate1D(all_tau_pts, all_T_pts, tY_loc, 2*n_tau_pts, 0, false, 2);
			double muY_loc = interpolate1D(all_tau_pts, all_mu_pts, tY_loc, 2*n_tau_pts, 0, false, 2);
			complex<double> factor_Y = (*Gtilde_Y)(-k, tY_loc) / tY_loc;

			double min_tp_tpp = min(tX_loc, tY_loc);
			double eta_at_min_tp_tpp = interpolate1D(all_tau_pts, running_integral_array, min_tp_tpp, 2*n_tau_pts, 0, false, 2);

			double sum_XY = exp(-abs(tX_loc - tY_loc) / tauC) * eta_at_min_tp_tpp / (2.0*tauC);
			sum_X += hw_loc * x_wts[ix] * factor_Y * sum_XY;
		}
		locsum += all_tau_wts[itp] * factor_X * sum_X;
	}

	delete [] x_pts;
	delete [] x_wts;

	return (locsum);
}

inline complex<double> Ctilde_n_n(double k, bool use_lattice_transport, double two_pi_DQ_T)
{
	complex<double> sum(0,0);
	sum = tau_integration(Gtilde_n, Gtilde_n, k, use_lattice_transport, two_pi_DQ_T);

	return ( sum );	//fm^2
}

////////////////////////////////////////////////////////////////////////////////
// Functions/stuff to solve for time-dependence of T and/or Tf
////////////////////////////////////////////////////////////////////////////////

struct rparams
{
	double tau;
};

//compute final time-step Tf and muf

int input_get_Tf_f (const gsl_vector * x, void * params, gsl_vector * f)
{
	const double x0 = gsl_vector_get (x, 0);	//T

	const double y0 = P(x0);						//defines P==0 curve

	gsl_vector_set (f, 0, y0);

	return GSL_SUCCESS;
}

void compute_Tf()
{
	const gsl_multiroot_fsolver_type *gsl_T;
	gsl_multiroot_fsolver *gsl_s;

	int status;
	size_t iter = 0;

	const size_t n = 1;
	struct rparams p = {taui};
	gsl_multiroot_function f = {&input_get_Tf_f, n, &p};

	double x_init[n] = {Ti};
	gsl_vector *x = gsl_vector_alloc (n);

	gsl_vector_set (x, 0, x_init[0]);

	gsl_T = gsl_multiroot_fsolver_hybrids;
	gsl_s = gsl_multiroot_fsolver_alloc (gsl_T, n);
	gsl_multiroot_fsolver_set (gsl_s, &f, x);

	do
	{
		iter++;

		status = gsl_multiroot_fsolver_iterate (gsl_s);

		if (status)   /* check if solver is stuck */
			break;

		status = gsl_multiroot_test_residual (gsl_s->f, 1e-7);
	}
	while (status == GSL_CONTINUE && iter < 1000);

	//finally, store results
	Tf = gsl_vector_get (gsl_s->x, 0);

	gsl_multiroot_fsolver_free (gsl_s);
	gsl_vector_free (x);

	return;
}

//compute T at each time step

int input_f (const gsl_vector * x, void * params, gsl_vector * f)
{
	double tau_local = ((struct rparams *) params)->tau;

	const double x0 = gsl_vector_get (x, 0);	//T

	const double y0 = s_vs_T(x0) - s_vs_tau(tau_local);

	gsl_vector_set (f, 0, y0);

	return GSL_SUCCESS;
}

void populate_T_vs_tau()
{
	for (int it = 0; it < n_tau_pts; ++it)
	{
		const gsl_multiroot_fsolver_type *gsl_T;
		gsl_multiroot_fsolver *gsl_s;

		int status;
		size_t iter = 0;

		const size_t n = 1;
		struct rparams p = {tau_pts[it]};
		gsl_multiroot_function f = {&input_f, n, &p};

		double x_init[n] = {guess_T(tau_pts[it])};
		gsl_vector *x = gsl_vector_alloc (n);

		gsl_vector_set (x, 0, x_init[0]);

		gsl_T = gsl_multiroot_fsolver_hybrids;
		gsl_s = gsl_multiroot_fsolver_alloc (gsl_T, n);
		gsl_multiroot_fsolver_set (gsl_s, &f, x);

		do
		{
			iter++;
			status = gsl_multiroot_fsolver_iterate (gsl_s);

			if (status)   /* check if solver is stuck */
				break;

			status = gsl_multiroot_test_residual (gsl_s->f, 1e-7);
		}
		while (status == GSL_CONTINUE && iter < 1000);

		//finally, store results
		T_pts[it] = gsl_vector_get (gsl_s->x, 0);

		gsl_multiroot_fsolver_free (gsl_s);
		gsl_vector_free (x);
	}
	
	return;
}

double compute_newTf(double new_tauf)
{
	const gsl_multiroot_fsolver_type *gsl_T;
	gsl_multiroot_fsolver *gsl_s;

	int status;
	size_t iter = 0;

	const size_t n = 1;
	struct rparams p = {new_tauf};
	gsl_multiroot_function f = {&input_f, n, &p};

	double x_init[n] = {guess_T(new_tauf)};
	gsl_vector *x = gsl_vector_alloc (n);

	gsl_vector_set (x, 0, x_init[0]);

	gsl_T = gsl_multiroot_fsolver_hybrids;
	gsl_s = gsl_multiroot_fsolver_alloc (gsl_T, n);
	gsl_multiroot_fsolver_set (gsl_s, &f, x);

	do
	{
		iter++;
		status = gsl_multiroot_fsolver_iterate (gsl_s);

		if (status)   /* check if solver is stuck */
			break;

		status = gsl_multiroot_test_residual (gsl_s->f, 1e-7);
	}
	while (status == GSL_CONTINUE && iter < 1000);

	//finally, store result
	double result = gsl_vector_get (gsl_s->x, 0);

	gsl_multiroot_fsolver_free (gsl_s);
	gsl_vector_free (x);
	
	return (result);
}

// End of file

#endif
