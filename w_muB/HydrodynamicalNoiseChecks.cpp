#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <complex>

//#include <gsl/gsl_integration.h>

using namespace std;

#include "defs.h"
#include "gauss_quadrature.h"

bool do_1p_calc = false;
bool do_HBT_calc = true;
bool scale_out_y_dependence = true;

const double infinity = 15.0;
double vs, Neff, tauf, taui, Tf, Ti, nu, nuVB, ds, A, m, sf;
double mByT, alpha0, phi0;
double chi_tilde_mu_mu, chi_tilde_T_mu, chi_tilde_T_T, Delta;

double exp_delta, exp_gamma, exp_nu;
double T0, mu0, Tc, Pc, nc, sc, wc, muc;
double A0, A2, A4, C0, B, mui, muf, xi0, xibar0, etaBYs, RD, sPERn, Nf, qD;

const int n_xi_pts = 200;
const int n_k_pts = 200;
const int n_tau_pts = 200;
double * xi_pts_0_inf, * xi_wts_0_inf, * xi_pts_minf_inf, * xi_wts_minf_inf;
double * k_pts, * k_wts;
double * tau_pts, * tau_wts;

int main(int argc, char *argv[])
{
	//set parameters corresponding to all different possible trajectories
	//everything in powers of fm
	double ctmms[3] = {5431.79, 5532.66, 5652.37};
	double ctTms[3] = {12739., 17872.4, 22182.7};
	double ctTTs[3] = {1.00076e6, 970752., 934775.};
	double sPERns[3] = {37.9844, 26.0847, 20.0631};
	double taufs[3] = {3.04, 3.30, 3.68};
	double sis[3] = {29.7527, 31.2566, 33.3389};	//fm^3
	double nis[3] = {0.783289, 1.19827, 1.6617};	//fm^3
	double Tfs[3] = {152.478, 149.887, 146.722};
	double mufs[3] = {187.98, 268.288, 340.174};
	double muis[3] = {420.0, 620.0, 820.0};

	//int chosen_trajectory = 1;	//numbering convention in paper
	int chosen_trajectory = atoi(argv[1]);
	chosen_trajectory--;		//use for array indexing throughout

    // initialize parameters
	Nf = 2.0;
	T0 = 170.0;
	mu0 = 1218.49;
	A4 = M_PI*M_PI*(16.0 + 10.5*Nf) / 90.0;
	A2 = Nf / 18.0;
	A0 = Nf / (324.0 * M_PI * M_PI);
	C0 = mu0*mu0*( A2 - 2.0*A0*mu0*mu0 / (T0*T0) );
	B = 0.8 * pow(T0, 4.0);
	Ti = 250.0;
	Tf = Tfs[chosen_trajectory];
	mui = muis[chosen_trajectory];
	muf = mufs[chosen_trajectory];
	taui = 0.5;		//fm/c
	tauf = taufs[chosen_trajectory];
    nu = 1.0 / (3.0 * M_PI);
	nuVB = nu;
    ds = 2.0;
    A = 1.0;
    m = 139.57;
	xi0 = 0.37;
	xibar0 = 0.69;
	qD = 1.0/xi0;
	RD = 1.05;
	mByT = m / Tf;
	sPERn = sPERns[chosen_trajectory];
	sf = s(Tf, muf);
	exp_delta = 4.815;
	exp_gamma = 1.24;
	exp_nu = 0.63;
	etaBYs = 1.0 / (4.0 * M_PI);
	si = sis[chosen_trajectory];
	ni = nis[chosen_trajectory];

	//set critical point quantities
	Tc = 160.0;
	muc = 411.74;
	Pc = P(Tc, muc);
	sc = s(Tc, muc);
	nc = n(Tc, muc);
	wc = w(Tc, muc);
	
	//set the susceptibilities
	double chi_mu_mu = ctmms[chosen_trajectory];
	double chi_T_mu = ctTms[chosen_trajectory];
	double chi_T_T = ctTTs[chosen_trajectory];
	Delta = chi_mu_mu * chi_T_T - chi_T_mu * chi_T_mu;
	chi_tilde_mu_mu = chi_mu_mu / Delta;
	chi_tilde_T_mu = chi_T_mu / Delta;
	chi_tilde_T_T = chi_T_T / Delta;

	double inv_vs2 = 1.0 / (vs*vs);

	double kappa = ds*A*tauf*Tf*Tf*Tf / (4.0*M_PI*M_PI);
	double kappap = ( (45.0 * ds * nu) / (4.0 * pow(M_PI, 4.0)*Neff*Tf*tauf) ) * pow( Ti/Tf, 2.0*(inv_vs2 - 2.0) );

    // set up grid points for integrations
    xi_pts_0_inf = new double [n_xi_pts];
    xi_pts_minf_inf = new double [n_xi_pts];
    tau_pts = new double [n_tau_pts];
    k_pts = new double [n_k_pts];
    xi_wts_0_inf = new double [n_xi_pts];
    xi_wts_minf_inf = new double [n_xi_pts];
    tau_wts = new double [n_tau_pts];
    k_wts = new double [n_k_pts];

    int tmp = gauss_quadrature(n_xi_pts, 1, 0.0, 0.0, 0.0, infinity, xi_pts_0_inf, xi_wts_0_inf);
    tmp = gauss_quadrature(n_xi_pts, 1, 0.0, 0.0, -infinity, infinity, xi_pts_minf_inf, xi_wts_minf_inf);
    tmp = gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, taui, tauf, tau_pts, tau_wts);
    tmp = gauss_quadrature(n_k_pts, 1, 0.0, 0.0, -infinity, infinity, k_pts, k_wts);

	//check Delta_lambda
	for (int it = 1; it < 41; ++it)
	{
		double t_loc = (double)it * 0.1;
		double T_loc = T(t_loc);
		double mu_loc = mu(t_loc);
		cout << chosen_trajectory + 1 << "   " << T_loc << "   " << mu_loc << "   " << t_loc << "   " << (1e-6)*Delta_lambda(T_loc, mu_loc) << endl;
	}

	if (1) return (0);

	if (do_1p_calc)
	{
		double norm = integrate_1D(norm_int, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts);
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
				sum += k_wts[ik] * exp(i * k * Delta_y)
						* ( Fts * conj(Fts) * Ctilde_s_s(k)
							+ Fts * conj(Fto) * Ctilde_s_omega(k)
							+ Fts * conj(Ftn) * Ctilde_s_n(k)
							+ Fto * conj(Fts) * Ctilde_omega_s(k)
							+ Fto * conj(Fto) * Ctilde_omega_omega(k)
							+ Fto * conj(Ftn) * Ctilde_omega_n(k)
							+ Ftn * conj(Fts) * Ctilde_n_s(k)
							+ Ftn * conj(Fto) * Ctilde_n_omega(k)
							+ Ftn * conj(Ftn) * Ctilde_n_n(k) );
			}

			complex<double> result = (kappa / (2.0 *M_PI * kappap * norm)) * sum;
			cout << Delta_y << "   " << result.real() << endl;
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
