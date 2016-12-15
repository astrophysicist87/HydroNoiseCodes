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

bool do_1p_calc;
bool do_HBT_calc;
bool scale_out_y_dependence = true;

bool white_noise, maxwell_cattaneo_noise;
int noise_mode = 0;		// 0 - white noise
						// 1 - Maxwell-Cattaneo noise

const double hbarC = 197.33;
const double k_infinity = 10.0;
const double xi_infinity = 4.0;
double vs, Neff, tauf, tau0, Tf, T0, nu, nuVB, ds, m;
double mByT, alpha, alpha0, phi0;
double CN_tau_D;

const int n_xi_pts = 100;
const int n_k_pts = 100;
const int n_tau_pts = 500;
double * xi_pts_0_inf, * xi_wts_0_inf, * xi_pts_minf_inf, * xi_wts_minf_inf;
double * k_pts, * k_wts;
double * tau_pts, * tau_wts;

double * tau_integral_factor_pts;

int main(int argc, char *argv[])
{
    // initialize parameters
    vs = 1./sqrt(3.0);
	alpha = 0.5*(1.0 - vs*vs);
    Neff = 47.5;
    tauf = 10.0;
    Tf = 150.0 / hbarC;
    //Tf = 150.0;
    T0 = 600.0 / hbarC;
    //T0 = 600.0;
    nu = 1.0 / (3.0 * M_PI);
	nuVB = nu;
	//nuVB = 0.0;
    ds = 2.0;
    m = 139.57 / hbarC;
    //m = 139.57;
    tau0 = invTemperature(Tf);
	mByT = m / Tf;

	double inv_vs2 = 1.0 / (vs*vs);

	//double kappa = ds*tauf*Tf*Tf*Tf / (4.0*M_PI*M_PI*hbarC*hbarC*hbarC);	//fm^{-2}
	//double kappap = ( (45.0 * ds * nu*hbarC) / (4.0 * pow(M_PI, 4.0)*Neff*Tf*tauf) ) * pow( T0/Tf, 2.0*(inv_vs2 - 2.0) );	//dimensionless
	double kappa = ds*tauf*Tf*Tf*Tf / (4.0*M_PI*M_PI);	//fm^{-2}
	double kappap = ( (45.0 * ds * nu) / (4.0 * pow(M_PI, 4.0)*Neff*Tf*tauf) ) * pow( T0/Tf, 2.0*(inv_vs2 - 2.0) );	//dimensionless

	do_1p_calc = (bool)atoi(argv[1]);
	do_HBT_calc = (bool)atoi(argv[2]);

	//set the noise mode and related quantities
	CN_tau_D = atof(argv[3]);	//fm
	//white_noise = bool(noise_mode == 0);
	//maxwell_cattaneo_noise = bool(noise_mode == 1);
	white_noise = bool(CN_tau_D < 1.e-6);
	//maxwell_cattaneo_noise = bool(CN_tau_D >= 1.e-6);
	maxwell_cattaneo_noise = not white_noise;

    // set up grid points for integrations
    xi_pts_0_inf = new double [n_xi_pts];
    xi_pts_minf_inf = new double [n_xi_pts];
    tau_pts = new double [n_tau_pts];
    k_pts = new double [n_k_pts];
    xi_wts_0_inf = new double [n_xi_pts];
    xi_wts_minf_inf = new double [n_xi_pts];
    tau_wts = new double [n_tau_pts];
    k_wts = new double [n_k_pts];

	tau_integral_factor_pts = new double [n_tau_pts];

    int tmp = gauss_quadrature(n_xi_pts, 1, 0.0, 0.0, 0.0, xi_infinity, xi_pts_0_inf, xi_wts_0_inf);
    tmp = gauss_quadrature(n_xi_pts, 1, 0.0, 0.0, -xi_infinity, xi_infinity, xi_pts_minf_inf, xi_wts_minf_inf);
    tmp = gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, tau0, tauf, tau_pts, tau_wts);
    tmp = gauss_quadrature(n_k_pts, 1, 0.0, 0.0, -k_infinity, k_infinity, k_pts, k_wts);

	set_tau_integral_factor_pts();

	if (do_1p_calc)
	{
		double norm = integrate_1D(norm_int, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts);
		vector<complex<double> > Ftr_vec, Fto_vec;
		vector<complex<double> > Ctrr_vec, Ctro_vec, Ctor_vec, Ctoo_vec;

		for (int ik = 0; ik < n_k_pts; ++ik)
		{
			double k = k_pts[ik];
			//current_ik = ik;
			Ftr_vec.push_back(Ftilde_rho(k));
			Fto_vec.push_back(Ftilde_omega(k));
			Ctrr_vec.push_back(Ctilde_rho_rho(k));
			Ctro_vec.push_back(Ctilde_rho_omega(k));
			Ctor_vec.push_back(Ctilde_omega_rho(k));
			Ctoo_vec.push_back(Ctilde_omega_omega(k));
		}

		for (int iDy = 0; iDy < 31; iDy++)
		{
			double Delta_y = (double)iDy * 0.2;
			complex<double> sum(0,0);
			for (int ik = 0; ik < n_k_pts; ++ik)
			{
				double k = k_pts[ik];
				//complex<double> Ftr = Ftilde_rho(k);
				//complex<double> Fto = Ftilde_omega(k);
				complex<double> Ftr = Ftr_vec[ik];
				complex<double> Fto = Fto_vec[ik];
				sum += k_wts[ik] * exp(i * k * Delta_y)
						* ( Ftr * conj(Ftr) * Ctrr_vec[ik]
							+ Ftr * conj(Fto) * Ctro_vec[ik]
							+ Fto * conj(Ftr) * Ctor_vec[ik]
							+ Fto * conj(Fto) * Ctoo_vec[ik] );
			}

			complex<double> result = (kappa / (2.0 *M_PI * kappap * norm)) * sum;
			cout << Delta_y << "   " << result.real() << endl;
		}
	}

	if (do_HBT_calc)
	{
		alpha0 = 1.0 / integrate_2D(alpha0_int, xi_pts_minf_inf, xi_pts_minf_inf, xi_wts_minf_inf, xi_wts_minf_inf, n_xi_pts, n_xi_pts);
		phi0 = integrate_2D(phi0_int, xi_pts_minf_inf, xi_pts_minf_inf, xi_wts_minf_inf, xi_wts_minf_inf, n_xi_pts, n_xi_pts);
		double norm = integrate_1D(norm_int, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts);
		for (int iDy = 0; iDy < 1; iDy++)
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
			double mean_R2l_vs_Dy = 0.5*tauf*tauf*phi0 / (cDy*cDy);
			double mean_R2l_vs_y1 = 0.5*tauf*tauf*phi0 / (cy1*cy1);
			double mean_R2l_vs_y2 = 0.5*tauf*tauf*phi0 / (cy2*cy2);
			complex<double> result = (ds*tauf*Tf*Tf*Tf / (4.0*M_PI*M_PI))*norm*scale_out_y_dep_factor * sum / (2.*M_PI*mean_R2l_vs_Dy);
			complex<double> result2 = (ds*tauf*Tf*Tf*Tf / (4.0*M_PI*M_PI))*norm*scale_out_y_dep_factor * sum / (2.*M_PI*mean_R2l_vs_y1*mean_R2l_vs_y2);
			cout << Delta_y << "   " << mean_R2l_vs_Dy << "   " << mean_R2l_vs_y1 << "   " << mean_R2l_vs_y2 << "   " << result.real() << "   " << result2.real() << endl;
		}
	}

	return 0;
}
