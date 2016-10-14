#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

using namespace std;

#include "defs.h"
#include "gauss_quadrature.h"

bool do_1p_calc = true;
bool do_HBT_calc = true;
const int particle_to_study = 1;	//1 is pion, 2 is proton

const double hbarC = 197.33;
const double k_infinity = 10.0;
const double xi_infinity = 4.0;
double vs, Neff, tauf, taui, Tf, Ti, nu, nuVB, ds, A, m, sf;
double mByT, alpha0, phi0;
double chi_tilde_mu_mu, chi_tilde_T_mu, chi_tilde_T_T, Delta;

double exp_delta, exp_gamma, exp_nu;
double T0, mu0, Tc, Pc, nc, sc, wc, muc;
double A0, A2, A4, C0, B, mui, muf, xi0, xibar0, etaBYs, RD, sPERn, Nf, qD, si, ni;
double a_at_tauf, vs2_at_tauf, vn2_at_tauf, vsigma2_at_tauf;
//complex<double> b_at_tauf;

const int n_xi_pts = 100;
const int n_k_pts = 100;
const int n_tau_pts = 100;
//int n_tau_pts;
double * xi_pts_minf_inf, * xi_wts_minf_inf;
double * k_pts, * k_wts;
double * tau_pts, * tau_wts;
double * T_pts, * mu_pts;

int main(int argc, char *argv[])
{
	//set parameters corresponding to all different possible trajectories
	Ti = 250.0;		//initial trajectory temperature
	double muis[3] = {420.0, 620.0, 820.0};

	//int chosen_trajectory = 1;	//numbering convention in paper
	//int chosen_trajectory = atoi(argv[1]);
	int chosen_trajectory = 3;
	int nt = atoi(argv[1]);
	chosen_trajectory--;		//use for array indexing throughout

	mui = muis[chosen_trajectory];

	set_phase_diagram_and_EOS_parameters();

	//other constants
	m = (particle_to_study == 1) ? 139.57 : 939.0;
	taui = 0.5;		//fm/c
	//taui = 0.493325;

	si = s(Ti, mui);
	ni = n(Ti, mui);
	sPERn = si / ni;

	compute_Tf_and_muf();
	//if (particle_to_study == 1) muf = 1e-6;

	sf = s(Tf, muf);
	tauf = si * taui / sf;

	set_critical_point_parameters();

    // initialize other parameters
    nu = 1.0 / (3.0 * M_PI);
	nuVB = nu;
    ds = 2.0;
    A = 1.0;
	mByT = m / Tf;

	//set the susceptibilities
	double chi_mu_mu = chi_mumu(Tf, muf);
	double chi_T_mu = chi_Tmu(Tf, muf);
	double chi_T_T = chi_TT(Tf, muf);
	Delta = chi_mu_mu * chi_T_T - chi_T_mu * chi_T_mu;
	chi_tilde_mu_mu = chi_mu_mu / Delta;
	chi_tilde_T_mu = chi_T_mu / Delta;
	chi_tilde_T_T = chi_T_T / Delta;

	//set parameters for FTd-Green's functions
	a_at_tauf = alpha(Tf, muf);
	vs2_at_tauf = vs2(Tf, muf);
	vn2_at_tauf = vn2(Tf, muf);
	vsigma2_at_tauf = vsigma2(Tf, muf);

    // set up grid points for integrations
    xi_pts_minf_inf = new double [n_xi_pts];
    tau_pts = new double [n_tau_pts];
    k_pts = new double [n_k_pts];
    xi_wts_minf_inf = new double [n_xi_pts];
    tau_wts = new double [n_tau_pts];
    k_wts = new double [n_k_pts];

    int tmp = gauss_quadrature(n_xi_pts, 1, 0.0, 0.0, -xi_infinity, xi_infinity, xi_pts_minf_inf, xi_wts_minf_inf);
    tmp = gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, taui, tauf, tau_pts, tau_wts);
    tmp = gauss_quadrature(n_k_pts, 1, 0.0, 0.0, -k_infinity, k_infinity, k_pts, k_wts);

	T_pts = new double [n_tau_pts];
	mu_pts = new double [n_tau_pts];

	//computes tau-dependence of T and mu for remainder of calculation
	populate_T_and_mu_vs_tau();

double Delta_t = (2.57 - 2.53) / (double)nt;

	//check Delta_lambda
	for (int it = 1; it < nt; ++it)
	{
		//double t_loc = tau_pts[it];
		//double T_loc = T_pts[it];
		//double mu_loc = mu_pts[it];
		double t_loc = 2.53 + (double)it * Delta_t;
		double T_loc = interpolate1D(tau_pts, T_pts, t_loc, n_tau_pts, 1, false);
		double mu_loc = interpolate1D(tau_pts, mu_pts, t_loc, n_tau_pts, 1, false);
		cout << setprecision(15) << t_loc << "   " << T_loc << "   " << mu_loc << "   " << (1e-6)*Delta_lambda(T_loc, mu_loc) << endl;
	}

	if (1) return (0);

	/*for (int ix = 0; ix < 200; ix++)
	{
		double x = -10.0 + (double)ix * 0.1;
		cout << x << "   " << Fs(x) << "   " << Fomega(x) << "   " << Fn(x) << endl;
	}
	if (1) return (0);*/

	if (do_1p_calc)
	{
		double norm = integrate_1D(norm_int, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts);
		//cout << "norm = " << setprecision(15) << norm << endl;
		for (int iDy = 0; iDy < 1; iDy++)
		{
			double Delta_y = (double)iDy * 0.1;
			complex<double> sum(0,0);
			for (int ik = 0; ik < n_k_pts; ++ik)
			{
				double k = k_pts[ik];
				complex<double> Fts = Ftilde_s(k);
				complex<double> Fto = Ftilde_omega(k);
				complex<double> Ftn = Ftilde_n(k);
				/*sum += k_wts[ik] * exp(i * k * Delta_y)
						* ( Fts * conj(Fts) * Ctilde_s_s(k)
							+ Fts * conj(Fto) * Ctilde_s_omega(k)
							+ Fts * conj(Ftn) * Ctilde_s_n(k)
							+ Fto * conj(Fts) * Ctilde_omega_s(k)
							+ Fto * conj(Fto) * Ctilde_omega_omega(k)
							+ Fto * conj(Ftn) * Ctilde_omega_n(k)
							+ Ftn * conj(Fts) * Ctilde_n_s(k)
							+ Ftn * conj(Fto) * Ctilde_n_omega(k)
							+ Ftn * conj(Ftn) * Ctilde_n_n(k) );*/
				complex<double> Ctss = Ctilde_s_s(k);
				complex<double> Ctso = Ctilde_s_omega(k);
				complex<double> Ctsn = Ctilde_s_n(k);
				complex<double> Ctos = Ctilde_omega_s(k);
				complex<double> Ctoo = Ctilde_omega_omega(k);
				complex<double> Cton = Ctilde_omega_n(k);
				complex<double> Ctns = Ctilde_n_s(k);
				complex<double> Ctno = Ctilde_n_omega(k);
				complex<double> Ctnn = Ctilde_n_n(k);
				sum += k_wts[ik] * exp(i * k * Delta_y)
						* ( Fts * conj(Fts) * Ctss
							+ Fts * conj(Fto) * Ctso
							+ Fts * conj(Ftn) * Ctsn
							+ Fto * conj(Fts) * Ctos
							+ Fto * conj(Fto) * Ctoo
							+ Fto * conj(Ftn) * Cton
							+ Ftn * conj(Fts) * Ctns
							+ Ftn * conj(Fto) * Ctno
							+ Ftn * conj(Ftn) * Ctnn );

				complex<double> contrib = k_wts[ik] * exp(i * k * Delta_y)
						* ( Fts * conj(Fts) * Ctss
							+ Fts * conj(Fto) * Ctso
							+ Fts * conj(Ftn) * Ctsn
							+ Fto * conj(Fts) * Ctos
							+ Fto * conj(Fto) * Ctoo
							+ Fto * conj(Ftn) * Cton
							+ Ftn * conj(Fts) * Ctns
							+ Ftn * conj(Fto) * Ctno
							+ Ftn * conj(Ftn) * Ctnn );

				/*cout << k_pts[ik] << "   " << Fts.real() << "   " << Fto.real() << "   " << Ftn.real() << "   "
						<< Ctss.real() << "   " << Ctso.real() << "   " << Ctsn.real() << "   "
						<< Ctos.real() << "   " << Ctoo.real() << "   " << Cton.real() << "   "
						<< Ctns.real() << "   " << Ctno.real() << "   " << Ctnn.real() << "   "
						<< Fts.imag() << "   " << Fto.imag() << "   " << Ftn.imag() << "   "
						<< Ctss.imag() << "   " << Ctso.imag() << "   " << Ctsn.imag() << "   "
						<< Ctos.imag() << "   " << Ctoo.imag() << "   " << Cton.imag() << "   "
						<< Ctns.imag() << "   " << Ctno.imag() << "   " << Ctnn.imag() << endl;*/
				//cout << k_pts[ik] << "   " << contrib.real() << "   " << sum.real() << endl;
			}

			sum *= hbarC*hbarC;  //hbarc^2 / hbarc

			complex<double> result = (exp(muf/Tf)*ds*tauf*Tf / hbarC / (2.0*M_PI*M_PI * norm)) * sum;	//this is the one I find
			//cout << n_tau_pts << "   " << Delta_y << "   " << result.real() /*<< "   " << sum.real() << "   " << exp(muf/Tf) << "   " << (tauf*Tf / hbarC) * ds*sum.real()/ (2.0*M_PI*M_PI * norm)*/ << endl;
		}
	}

	/*if (do_HBT_calc)
	{
		alpha0 = 1.0 / integrate_2D(alpha0_int, xi_pts_minf_inf, xi_pts_minf_inf, xi_wts_minf_inf, xi_wts_minf_inf, n_xi_pts, n_xi_pts);
		phi0 = integrate_2D(phi0_int, xi_pts_minf_inf, xi_pts_minf_inf, xi_wts_minf_inf, xi_wts_minf_inf, n_xi_pts, n_xi_pts);
		for (int iDy = 0; iDy < 101; iDy++)
		{
			double Delta_y = (double)iDy * 0.1;
			//option #1
			double y1 = Delta_y;
			double y2 = 0.0;
			//option #2
			//double y1 = 0.5*Delta_y;
			//double y2 = -0.5*Delta_y;
			//option #3
			//double y1 = Delta_y;
			//double y2 = Delta_y;
//cout << alpha0 << "   " << phi0 << endl;
			double cy1 = cosh(y1);
			double cy2 = cosh(y2);
			double cDy = cosh(y1-y2);
			double scale_out_y_dep_factor = 1.0;
			if (scale_out_y_dependence)
				scale_out_y_dep_factor = cy1*cy1*cy2*cy2 / (cDy*cDy);
			complex<double> sum(0,0);
			for (int ik = 0; ik < n_k_pts; ++ik)
			{
				double k = k_pts[ik];
				complex<double> gt_r = gt_rho(k);
				complex<double> gt_o = gt_omega(k);
				sum += k_wts[ik] * exp(i * k * Delta_y)
						* ( gt_r * conj(gt_r) * Ctilde_rho_rho(k)
							+ gt_r * conj(gt_o) * Ctilde_rho_omega(k)
							+ gt_o * conj(gt_r) * Ctilde_omega_rho(k)
							+ gt_o * conj(gt_o) * Ctilde_omega_omega(k) )
						/ (cy1*cy1*cy2*cy2);
//cout << k << "   " << gt_r * conj(gt_r) << "   " << gt_r * conj(gt_o) << "   " << gt_o * conj(gt_r) << "   " << gt_o * conj(gt_o) << "   "
//		<< Ctilde_rho_rho(k) << "   " << Ctilde_rho_omega(k)<< "   " << Ctilde_omega_rho(k) << "   " << Ctilde_omega_omega(k) << endl;
			}

			//complex<double> result = sum;
			double mean_R2l_vs_y = 0.5*tauf*tauf*phi0 / (cDy*cDy);
			complex<double> result = scale_out_y_dep_factor * sum / (2.*M_PI*mean_R2l_vs_y);
			cout << Delta_y << "   " << result.real() << endl;
		}
	}*/

	return 0;
}
