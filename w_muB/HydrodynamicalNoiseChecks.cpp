#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>
#include <algorithm>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

using namespace std;

#include "defs.h"
#include "gauss_quadrature.h"

bool do_1p_calc = true;
bool do_HBT_calc = true;
bool scale_out_y_dependence = false;
const int particle_to_study = 1;	//1 is pion, 2 is proton

const double hbarC = 197.33;
const double k_infinity = 10.0;
const double xi_infinity = 4.0;
double vs, Neff, tauf, taui, Tf, Ti, nu, nuVB, ds, m, sf;
double mByT, alpha0, psi0;
double chi_tilde_mu_mu, chi_tilde_T_mu, chi_tilde_T_T, Delta;

double exp_delta, exp_gamma, exp_nu;
double T0, mu0, Tc, Pc, nc, sc, wc, muc;
double A0, A2, A4, C0, B, mui, muf, xi0, xibar0, etaBYs, RD, sPERn, Nf, qD, si, ni;
double a_at_tauf, vs2_at_tauf, vn2_at_tauf, vsigma2_at_tauf;

const int n_xi_pts = 200;
const int n_k_pts = 200;
const int n_tau_pts = 200;
double * xi_pts_minf_inf, * xi_wts_minf_inf;
double * k_pts, * k_wts;
double * tau_pts, * tau_wts;
double * tau_pts_lower, * tau_wts_lower;
//double * tau_pts_middle, * tau_wts_middle;
double * tau_pts_upper, * tau_wts_upper;
double * T_pts, * mu_pts;
double * T_pts_lower, * mu_pts_lower;
//double * T_pts_middle, * mu_pts_middle;
double * T_pts_upper, * mu_pts_upper;

int main(int argc, char *argv[])
{
	//set parameters corresponding to all different possible trajectories
	//Ti = 250.0;		//initial trajectory temperature
	Ti = 250.0 / hbarC;		//initial trajectory temperature
	double muis[3] = {420.0, 620.0, 820.0};

	int chosen_trajectory = atoi(argv[1]);
	chosen_trajectory--;		//use for array indexing throughout

	mui = muis[chosen_trajectory];
	mui /= hbarC;

	set_phase_diagram_and_EOS_parameters();

	//other constants
	//m = (particle_to_study == 1) ? 139.57 : 939.0;
	m = (particle_to_study == 1) ? 139.57 / hbarC : 939.0 / hbarC;
	taui = 0.5;		//fm/c

	si = s(Ti, mui);
	ni = n(Ti, mui);
	sPERn = si / ni;

	compute_Tf_and_muf();
	//if (particle_to_study == 1) muf = 1e-6 / hbarC;

	sf = s(Tf, muf);
	tauf = si * taui / sf;

	set_critical_point_parameters();

	if (particle_to_study == 1) muf = 1e-6 / hbarC;

    // initialize other parameters
    nu = 1.0 / (3.0 * M_PI);
	nuVB = nu;
    ds = 2.0;
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
    tau_pts_lower = new double [n_tau_pts];
    tau_pts_upper = new double [n_tau_pts];
    tau_pts = new double [n_tau_pts];
    k_pts = new double [n_k_pts];
    xi_wts_minf_inf = new double [n_xi_pts];
    tau_wts_lower = new double [n_tau_pts];
    tau_wts_upper = new double [n_tau_pts];
    tau_wts = new double [n_tau_pts];
    k_wts = new double [n_k_pts];

    int tmp = gauss_quadrature(n_xi_pts, 1, 0.0, 0.0, -xi_infinity, xi_infinity, xi_pts_minf_inf, xi_wts_minf_inf);
    tmp = gauss_quadrature(n_k_pts, 1, 0.0, 0.0, -k_infinity, k_infinity, k_pts, k_wts);
    tmp = gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, taui, tauf, tau_pts, tau_wts);

	T_pts_lower = new double [n_tau_pts];
	T_pts_upper = new double [n_tau_pts];
	T_pts = new double [n_tau_pts];
	mu_pts_lower = new double [n_tau_pts];
	mu_pts_upper = new double [n_tau_pts];
	mu_pts = new double [n_tau_pts];

	//computes tau-dependence of T and mu for remainder of calculation
	populate_T_and_mu_vs_tau();

	int nt = 1000000;
	double max_DL, tauc, width;
	break_up_integral(nt, max_DL, tauc, width);
    tmp = gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, taui, tauc, tau_pts_lower, tau_wts_lower);
    tmp = gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, tauc, tauf, tau_pts_upper, tau_wts_upper);
	populate_T_and_mu_vs_tau_part2();

	//if (1) return (0);

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
			}

			//complex<double> result = (exp(muf/Tf)*ds*tauf*Tf / hbarC / (2.0*M_PI*M_PI * norm)) * sum;
			complex<double> result = (exp(muf/Tf)*ds*tauf*Tf / (2.0*M_PI*M_PI * norm)) * sum;
			cout << Delta_y << "   " << result.real() << endl;
		}
	}

	if (do_HBT_calc)
	{
		double norm = integrate_1D(norm_int, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts);
		alpha0 = get_alpha0();
		psi0 = get_psi0();
		for (int iDy = 0; iDy < 1; iDy++)
		{
			double Delta_y = (double)iDy * 0.1;
			//option #1
			double y1 = Delta_y;
			double y2 = 0.0;
			//option #2
			//double y1 = 0.5*Delta_y;
			//double y2 = -0.5*Delta_y;

			//scale out y-dependence, if desired
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

				complex<double> Fbts = Fbtilde_s(k);
				complex<double> Fbto = Fbtilde_omega(k);
				complex<double> Fbtn = Fbtilde_n(k);

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
						* ( Fbts * conj(Fbts) * Ctss
							+ Fbts * conj(Fbto) * Ctso
							+ Fbts * conj(Fbtn) * Ctsn
							+ Fbto * conj(Fbts) * Ctos
							+ Fbto * conj(Fbto) * Ctoo
							+ Fbto * conj(Fbtn) * Cton
							+ Fbtn * conj(Fbts) * Ctns
							+ Fbtn * conj(Fbto) * Ctno
							+ Fbtn * conj(Fbtn) * Ctnn );
			}

			sum *= pow(tauf, 4.0)/ (cy1*cy1*cy2*cy2);

			double mean_R2l_vs_Dy = 0.5*tauf*tauf*psi0 / (cDy*cDy);
			double mean_R2l_vs_y1 = 0.5*tauf*tauf*psi0 / (cy1*cy1);
			double mean_R2l_vs_y2 = 0.5*tauf*tauf*psi0 / (cy2*cy2);
			complex<double> result = (exp(muf/Tf)*ds*tauf*Tf*norm / (2.0*M_PI*M_PI))*scale_out_y_dep_factor * sum / (mean_R2l_vs_Dy);
			complex<double> result2 = (exp(muf/Tf)*ds*tauf*Tf*norm / (2.0*M_PI*M_PI))*scale_out_y_dep_factor * sum / (mean_R2l_vs_y1*mean_R2l_vs_y2);
			cout << Delta_y << "   " << mean_R2l_vs_Dy << "   " << mean_R2l_vs_y1 << "   " << mean_R2l_vs_y2 << "   " << result.real() << "   " << result2.real() << endl;
		}
	}

	return 0;
}
