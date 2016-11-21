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

//bool do_1p_calc = true;
//bool do_HBT_calc = true;
bool do_1p_calc;
bool do_HBT_calc;
bool scale_out_y_dependence = false;
bool white_noise, maxwell_cattaneo_noise, gurtin_pipkin_noise;
int noise_mode = 1;		// 0 - white noise
						// 1 - Maxwell-Cattaneo noise
						// 2 - Gurtin-Pipkin noise

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

double screw_it_up_factor;
int current_ik;
const double mu_pion = 1.e-3 / hbarC;
double mu_proton, mu_part;

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
	double screw_it_up_factors[3] = {0.0, 1.0, 2.0};
	current_ik = 0;

	int chosen_trajectory = atoi(argv[1]);
	chosen_trajectory--;		//use for array indexing throughout
	do_1p_calc = (bool)atoi(argv[2]);
	do_HBT_calc = (bool)atoi(argv[3]);

	//set the noise mode
	white_noise = bool(noise_mode == 0);
	maxwell_cattaneo_noise = bool(noise_mode == 1);
	gutrin_pipkin_noise = bool(noise_mode == 2);

	//use this to see why all three trajectories zero at roughly same Delta_y point
	screw_it_up_factor = screw_it_up_factors[chosen_trajectory];

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

	sf = s(Tf, muf);
	tauf = si * taui / sf;

	set_critical_point_parameters();

	mu_proton = muf;
	mu_part = (particle_to_study == 1) ? mu_pion : mu_proton;
	//mu_part = muf;

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

//cout << setprecision(15) << taui << "   " << tauf << endl;


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

	//for (int it  = 0; it < n_tau_pts; it++)
	//	cout << setprecision(15) << it << "   " << tau_pts[it]*1000.0/hbarC << "   " << T_pts[it]*hbarC / 1000.0 << "   " << mu_pts[it]*hbarC / 1000.0 << endl;

	//if (1) return (0);

	/*cout << hbarC*Tf << "   " << hbarC*muf << "   " << hbarC*hbarC*hbarC*sf << "   " << chi_tilde_mu_mu/(hbarC*hbarC) << "   " << chi_tilde_T_mu/(hbarC*hbarC) << "   " << chi_tilde_T_T/(hbarC*hbarC) << endl << endl;

	for (int ix = 0; ix < 21; ix++)
	{
		double x = -1.0 + (double)ix * 0.1;
		cout << x << "   " << Fs(x) << "   " << Fomega(x) << "   " << Fn(x) << endl;
	}
	if (1) return (0);*/

	if (do_1p_calc)
	{
		double norm = integrate_1D(norm_int, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts);
		//cout << "norm = " << setprecision(15) << norm << endl;

		vector<complex<double> > Fts_vec, Fto_vec, Ftn_vec;
		vector<complex<double> > Ctss_vec, Ctso_vec, Ctsn_vec, Ctos_vec, Ctoo_vec, Cton_vec, Ctns_vec, Ctno_vec, Ctnn_vec;

		if (gurtin_pipkin_noise)
			set_noise_integral_points();

		for (int ik = 0; ik < n_k_pts; ++ik)
		{
			double k = k_pts[ik];
			current_ik = ik;
			Fts_vec.push_back(Ftilde_s(k));
			Fto_vec.push_back(Ftilde_omega(k));
			Ftn_vec.push_back(Ftilde_n(k));
			Ctss_vec.push_back(Ctilde_s_s(k));
			Ctso_vec.push_back(Ctilde_s_omega(k));
			Ctsn_vec.push_back(Ctilde_s_n(k));
			Ctos_vec.push_back(Ctilde_omega_s(k));
			Ctoo_vec.push_back(Ctilde_omega_omega(k));
			Cton_vec.push_back(Ctilde_omega_n(k));
			Ctns_vec.push_back(Ctilde_n_s(k));
			Ctno_vec.push_back(Ctilde_n_omega(k));
			Ctnn_vec.push_back(Ctilde_n_n(k));
		}

		for (int iDy = 0; iDy < 51; iDy++)
		{
			double Delta_y = (double)iDy * 0.1;
			complex<double> sum(0,0);
			complex<double> sum00(0,0), sum01(0,0), sum02(0,0), sum10(0,0), sum11(0,0), sum12(0,0), sum20(0,0), sum21(0,0), sum22(0,0);
			for (int ik = 0; ik < n_k_pts; ++ik)
			{
				double k = k_pts[ik];
				/*complex<double> Fts = Ftilde_s(k);
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
				complex<double> Ctnn = Ctilde_n_n(k);*/
				complex<double> Fts = Fts_vec[ik];
				complex<double> Fto = Fto_vec[ik];
				complex<double> Ftn = Ftn_vec[ik];
				complex<double> Ctss = Ctss_vec[ik];
				complex<double> Ctso = Ctso_vec[ik];
				complex<double> Ctsn = Ctsn_vec[ik];
				complex<double> Ctos = Ctos_vec[ik];
				complex<double> Ctoo = Ctoo_vec[ik];
				complex<double> Cton = Cton_vec[ik];
				complex<double> Ctns = Ctns_vec[ik];
				complex<double> Ctno = Ctno_vec[ik];
				complex<double> Ctnn = Ctnn_vec[ik];

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
				sum00 += k_wts[ik] * exp(i * k * Delta_y) * Fts * conj(Fts) * Ctss;
				sum01 += k_wts[ik] * exp(i * k * Delta_y) * Fts * conj(Fto) * Ctso;
				sum02 += k_wts[ik] * exp(i * k * Delta_y) * Fts * conj(Ftn) * Ctsn;
				sum10 += k_wts[ik] * exp(i * k * Delta_y) * Fto * conj(Fts) * Ctos;
				sum11 += k_wts[ik] * exp(i * k * Delta_y) * Fto * conj(Fto) * Ctoo;
				sum12 += k_wts[ik] * exp(i * k * Delta_y) * Fto * conj(Ftn) * Cton;
				sum20 += k_wts[ik] * exp(i * k * Delta_y) * Ftn * conj(Fts) * Ctns;
				sum21 += k_wts[ik] * exp(i * k * Delta_y) * Ftn * conj(Fto) * Ctno;
				sum22 += k_wts[ik] * exp(i * k * Delta_y) * Ftn * conj(Ftn) * Ctnn;
			}

			//complex<double> result = (exp(muf/Tf)*ds*tauf*Tf / hbarC / (2.0*M_PI*M_PI * norm)) * sum;
			complex<double> result = (exp(mu_part/Tf)*ds*tauf*Tf / (2.0*M_PI*M_PI * norm)) * sum;
			//cerr << "SUMMARY: " << Delta_y << "   " << sum00.real() << "   " << sum01.real() << "   " << sum02.real() << "   "
			//		<< sum10.real() << "   " << sum11.real() << "   " << sum12.real() << "   " << sum20.real() << "   "
			//		<< sum21.real() << "   " << sum22.real() << "   " << sum.real() << "   " << result.real() << endl;
			cout << setprecision(15) << Delta_y << "   " << result.real() << endl;
		}
	}

	if (do_HBT_calc)
	{
		double norm = integrate_1D(norm_int, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts);
		alpha0 = get_alpha0();
		psi0 = get_psi0();

		vector<complex<double> > Fbts_vec, Fbto_vec, Fbtn_vec;
		vector<complex<double> > Ctss_vec, Ctso_vec, Ctsn_vec, Ctos_vec, Ctoo_vec, Cton_vec, Ctns_vec, Ctno_vec, Ctnn_vec;

		if (gurtin_pipkin_noise)
			set_noise_integral_points();

		for (int ik = 0; ik < n_k_pts; ++ik)
		{
			double k = k_pts[ik];
			current_ik = ik;
			Fbts_vec.push_back(Fbtilde_s(k));
			Fbto_vec.push_back(Fbtilde_omega(k));
			Fbtn_vec.push_back(Fbtilde_n(k));
			Ctss_vec.push_back(Ctilde_s_s(k));
			Ctso_vec.push_back(Ctilde_s_omega(k));
			Ctsn_vec.push_back(Ctilde_s_n(k));
			Ctos_vec.push_back(Ctilde_omega_s(k));
			Ctoo_vec.push_back(Ctilde_omega_omega(k));
			Cton_vec.push_back(Ctilde_omega_n(k));
			Ctns_vec.push_back(Ctilde_n_s(k));
			Ctno_vec.push_back(Ctilde_n_omega(k));
			Ctnn_vec.push_back(Ctilde_n_n(k));
		}

		for (int iDy = 0; iDy < 51; iDy++)
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

				complex<double> Fbts = Fbts_vec[ik];
				complex<double> Fbto = Fbto_vec[ik];
				complex<double> Fbtn = Fbtn_vec[ik];

				complex<double> Ctss = Ctss_vec[ik];
				complex<double> Ctso = Ctso_vec[ik];
				complex<double> Ctsn = Ctsn_vec[ik];
				complex<double> Ctos = Ctos_vec[ik];
				complex<double> Ctoo = Ctoo_vec[ik];
				complex<double> Cton = Cton_vec[ik];
				complex<double> Ctns = Ctns_vec[ik];
				complex<double> Ctno = Ctno_vec[ik];
				complex<double> Ctnn = Ctnn_vec[ik];

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

			double A = 1.0; //fm^2
			double dNdy = exp(mu_part/Tf)*ds*tauf*A*Tf*Tf*Tf*norm / (4.0*M_PI*M_PI);	//A used to calculate dN/dy, used nowhere else
			double mean_R2l_vs_Dy = 0.5*tauf*tauf*psi0 / (cDy*cDy);
			double mean_R2l_vs_y1 = 0.5*tauf*tauf*psi0 / (cy1*cy1);
			double mean_R2l_vs_y2 = 0.5*tauf*tauf*psi0 / (cy2*cy2);
			complex<double> result = (exp(mu_part/Tf)*ds*tauf*Tf*Tf*Tf*norm / (2.0*M_PI*M_PI))*scale_out_y_dep_factor * sum / (mean_R2l_vs_Dy);
			complex<double> result2 = (exp(mu_part/Tf)*ds*tauf*Tf*Tf*Tf*norm / (2.0*M_PI*M_PI))*scale_out_y_dep_factor * sum / (mean_R2l_vs_y1*mean_R2l_vs_y2);
			cout << Delta_y << "   " << dNdy << "   " << mean_R2l_vs_Dy << "   " << mean_R2l_vs_y1 << "   " << mean_R2l_vs_y2 << "   " << result.real() << "   " << result2.real() << endl;
		}
	}

	return 0;
}
