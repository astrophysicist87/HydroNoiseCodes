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

const double infinity = 10.0;
double vs, Neff, tauf, tau0, Tf, T0, nu, nuVB, ds, A, m;
double mByT, alpha, alpha0, phi0;

const int n_xi_pts = 80;
const int n_k_pts = 80;
const int n_tau_pts = 80;
double * xi_pts_0_inf, * xi_wts_0_inf, * xi_pts_minf_inf, * xi_wts_minf_inf;
double * k_pts, * k_wts;
double * tau_pts, * tau_wts;

int main()
{
    // initialize parameters
    vs = 1./sqrt(3.0);
	alpha = 0.5*(1.0 - vs*vs);
    Neff = 47.5;
    tauf = 10.0;
    Tf = 150.0 / 197.33;
    T0 = 600.0 / 197.33;
    nu = 1.0 / (3.0 * M_PI);
	nuVB = nu;
    ds = 2.0;
    A = 1.0;
    m = 139.57 / 197.33;
    tau0 = invTemperature(Tf);
	mByT = m / Tf;

	double inv_vs2 = 1.0 / (vs*vs);

	double kappa = ds*A*tauf*Tf*Tf*Tf / (4.0*M_PI*M_PI);
	double kappap = ( (45.0 * ds * nu) / (4.0 * pow(M_PI, 4.0)*Neff*Tf*tauf) ) * pow( T0/Tf, 2.0*(inv_vs2 - 2.0) );

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
    tmp = gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, tau0, tauf, tau_pts, tau_wts);
    tmp = gauss_quadrature(n_k_pts, 1, 0.0, 0.0, -infinity, infinity, k_pts, k_wts);

	if (do_1p_calc)
	{
		double norm = integrate_1D(norm_int, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts);
		for (int iDy = 0; iDy < 61; iDy++)
		{
			double Delta_y = (double)iDy * 0.1;
			complex<double> sum(0,0);
			for (int ik = 0; ik < n_k_pts; ++ik)
			{
				double k = k_pts[ik];
				complex<double> Ftr = Ftilde_rho(k);
				complex<double> Fto = Ftilde_omega(k);
				sum += k_wts[ik] * exp(i * k * Delta_y)
						* ( Ftr * conj(Ftr) * Ctilde_rho_rho(k)
							+ Ftr * conj(Fto) * Ctilde_rho_omega(k)
							+ Fto * conj(Ftr) * Ctilde_omega_rho(k)
							+ Fto * conj(Fto) * Ctilde_omega_omega(k) );
			}

			complex<double> result = (kappa / (2.0 *M_PI * kappap * norm)) * sum;
			//cout << Delta_y << "   " << result.real() << endl;
		}
	}

	if (do_HBT_calc)
	{
		alpha0 = 1.0 / integrate_2D(alpha_int, xi_pts_minf_inf, xi_pts_minf_inf, xi_wts_minf_inf, xi_wts_minf_inf, n_xi_pts, n_xi_pts);
		for (int iDy = 0; iDy < 61; iDy++)
		{
			double Delta_y = (double)iDy * 0.1;
			//double y1 = Delta_y;
			//double y2 = 0.0;
			double y1 = 0.5*Delta_y;
			double y2 = -0.5*Delta_y;
			complex<double> sum(0,0);
			for (int ik = 0; ik < n_k_pts; ++ik)
			{
				double k = k_pts[ik];
				complex<double> Ftr1 = Fbt_rho(y1, k);
				complex<double> Ftr2 = Fbt_rho(y2, k);
				complex<double> Fto1 = Fbt_omega(y1, k);
				complex<double> Fto2 = Fbt_omega(y2, k);
				//cout << k << "   " << Ftr1 << "   " << Ftr2 << "   " << Fto1 << "   " << Fto2 << endl;
				sum += k_wts[ik]
						* ( Ftr1 * conj(Ftr2) * Ctilde_rho_rho(k)
							+ Ftr1 * conj(Fto2) * Ctilde_rho_omega(k)
							+ Fto1 * conj(Ftr2) * Ctilde_omega_rho(k)
							+ Fto1 * conj(Fto2) * Ctilde_omega_omega(k) );
			}

			//complex<double> result = sum;
			double mean_R2l_vs_y = -0.25*tauf*tauf*get_phi0(2.0*y1);
			complex<double> result = (kappa / (2.0 *M_PI * kappap)) * sum / mean_R2l_vs_y;
			cout << Delta_y << "   " << sum.real() << "   " << mean_R2l_vs_y << "   " << result.real() << endl;
		}
	}

	return 0;
}
