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

//const int particle_to_study = 2
	// 1 - pion
	// 2 - proton
	// 3 - kaon
int particle_to_study;

const double hbarC = 197.33;
const double Cem = 2.0 / 3.0;	//my current best guess
const double k_infinity = 10.0;
const double xi_infinity = 4.0;
const int n_Dy = 51;

long n_interp;

double vs, Neff, tauf, taui, Tf, Ti, nu, nuVB, ds, m, sf;
double mByT, alpha0, psi0;
double chi_tilde_mu_mu, chi_tilde_T_mu, chi_tilde_T_T, chi_mu_mu, chi_T_mu, chi_T_T, Delta;

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
double * T_pts;

double * interp_T_pts, * interp_transport_pts;

int main(int argc, char *argv[])
{
	//set parameters corresponding to all different possible trajectories
	particle_to_study = atoi(argv[1]);
	Ti = atoi(argv[2]) / hbarC;		//initial trajectory temperature
	fraction_of_evolution = atof(argv[3]);	//what "snapshot" to take of evolution - 1.0 means the full evolution as usual

	set_phase_diagram_and_EOS_parameters();

	//other constants
	double Delta_y_step = 0.0;
	switch (particle_to_study)
	{
		case 1:	//pion
			m = 139.57 / hbarC;
			Delta_y_step = 0.1;
			break;
		case 2:	//proton
			m = 939.0 / hbarC;
			Delta_y_step = 0.05;
			break;
		case 3:	//kaon
			m = 493.68 / hbarC;
			Delta_y_step = 0.075;
			break;
		default:
			cerr << "Not a supported particle!" << endl;
			exit(1);
			break;
	}

	taui = 0.5;		//fm/c

	si = s_vs_T(Ti);

	compute_Tf();

	sf = s_vs_T(Tf);
	tauf = si * taui / sf;

	// added this to allow "snapshots" throughout evolution
	tauf = taui + fraction_of_evolution * (tauf - taui);
	sf = s_vs_tau(tauf);
	Tf = compute_newTf(tauf);
	// rest of code runs the same as before

    // initialize other parameters
    ds = 2.0;
	mByT = m / Tf;

	//set the susceptibilities
	chi_mu_mu = chi_mumu(Tf);
	chi_T_mu = chi_Tmu(Tf);
	chi_T_T = chi_TT(Tf);
	Delta = chi_mu_mu * chi_T_T - chi_T_mu * chi_T_mu;
	chi_tilde_mu_mu = chi_mu_mu / Delta;
	chi_tilde_T_mu = chi_T_mu / Delta;
	chi_tilde_T_T = chi_T_T / Delta;

	//output some parameters from calculation
	cout 	<< "#########################################################" << endl
			<< "# Using following parameters:" << endl
			<< "# taui = " << taui << " fm/c, tauf = " << tauf << " fm/c" << endl
			<< "# si = " << si << ", sf = " << sf << endl
			<< "# Ti = " << Ti*hbarC << " MeV, Tf = " << Tf*hbarC << " MeV" << endl
			<< "# chi_{T,T} = " << chi_T_T << ", chi_{T,mu} = chi_{mu,T} = " << chi_T_mu << ", chi_{mu,mu} = " << chi_mu_mu << endl
			<< "#########################################################" << endl;

    // set up grid points for integrations
    xi_pts_minf_inf = new double [n_xi_pts];
    tau_pts = new double [n_tau_pts];
    k_pts = new double [n_k_pts];
    xi_wts_minf_inf = new double [n_xi_pts];
    tau_wts = new double [n_tau_pts];
    k_wts = new double [n_k_pts];

    int tmp = gauss_quadrature(n_xi_pts, 1, 0.0, 0.0, -xi_infinity, xi_infinity, xi_pts_minf_inf, xi_wts_minf_inf);
    tmp = gauss_quadrature(n_k_pts, 1, 0.0, 0.0, -k_infinity, k_infinity, k_pts, k_wts);
    tmp = gauss_quadrature(n_tau_pts, 1, 0.0, 0.0, taui, tauf, tau_pts, tau_wts);

	T_pts = new double [n_tau_pts];

	//computes tau-dependence of T and mu for remainder of calculation
	populate_T_vs_tau();

	// read in table for interpolation of transport coefficients
	n_interp = get_filelength("transport.dat");
	interp_T_pts = new double [n_interp];
	interp_transport_pts = new double [n_interp];
	read_in_transport_table(n_interp, "transport.dat");

	/*for (int it = 0; it < n_tau_pts; it++)
	{
			double T_loc = T_pts[it];
			double sigma_By_Cem_T = interpolate1D(interp_T_pts, interp_transport_pts, T_loc, n_interp, 0, 1);
			double result = Cem * T_loc * sigma_By_Cem_T;
		cout << setprecision(15) << it << "   " << tau_pts[it] << "   " << T_pts[it] << "   " << 0.0 << "   " << result << endl;
	}
	if (1) return (0);*/

	double norm = integrate_1D(norm_int, xi_pts_minf_inf, xi_wts_minf_inf, n_xi_pts);

	vector<complex<double> > Ftn_vec;
	vector<complex<double> > Ctnn_lattice_vec;
	vector<complex<double> > Ctnn_2piDQT_05_vec;
	vector<complex<double> > Ctnn_2piDQT_10_vec;
	vector<complex<double> > Ctnn_2piDQT_15_vec;

	for (int ik = 0; ik < n_k_pts; ++ik)
	{
		double k = k_pts[ik];
		Ftn_vec.push_back(Ftilde_n(k));
		Ctnn_lattice_vec.push_back(Ctilde_n_n(k, true, 0.0));
		Ctnn_2piDQT_05_vec.push_back(Ctilde_n_n(k, false, 0.5));
		Ctnn_2piDQT_10_vec.push_back(Ctilde_n_n(k, false, 1.0));
		Ctnn_2piDQT_15_vec.push_back(Ctilde_n_n(k, false, 1.5));
	}

	complex<double> intercept_lattice(0,0), intercept_2piDQT_05(0,0), intercept_2piDQT_10(0,0), intercept_2piDQT_15(0,0);
	for (int iDy = 0; iDy < n_Dy; iDy++)
	{
		double Delta_y = (double)iDy * Delta_y_step;

		complex<double> sum_lattice(0,0), sum_2piDQT_05(0,0), sum_2piDQT_10(0,0), sum_2piDQT_15(0,0);
		for (int ik = 0; ik < n_k_pts; ++ik)
		{
			double k = k_pts[ik];

			complex<double> Ftn = Ftn_vec[ik];
			complex<double> Ctnn_lattice = Ctnn_lattice_vec[ik];
			complex<double> Ctnn_2piDQT_05 = Ctnn_2piDQT_05_vec[ik];
			complex<double> Ctnn_2piDQT_10 = Ctnn_2piDQT_10_vec[ik];
			complex<double> Ctnn_2piDQT_15 = Ctnn_2piDQT_15_vec[ik];

			sum_lattice += k_wts[ik] * exp(i * k * Delta_y)
					* ( Ftn * conj(Ftn) * Ctnn_lattice );
			sum_2piDQT_05 += k_wts[ik] * exp(i * k * Delta_y)
					* ( Ftn * conj(Ftn) * Ctnn_2piDQT_05 );
			sum_2piDQT_10 += k_wts[ik] * exp(i * k * Delta_y)
					* ( Ftn * conj(Ftn) * Ctnn_2piDQT_10 );
			sum_2piDQT_15 += k_wts[ik] * exp(i * k * Delta_y)
					* ( Ftn * conj(Ftn) * Ctnn_2piDQT_15 );
		}

		//adjust average value to get true charge balance function
		if (iDy == 0)	//assume we start at Delta y == 0
		{
			intercept_lattice = sum_lattice;
			intercept_2piDQT_05 = sum_2piDQT_05;
			intercept_2piDQT_10 = sum_2piDQT_10;
			intercept_2piDQT_15 = sum_2piDQT_15;
		}

		//N.B. - factors of 2 convert to charge balance function
		complex<double> result_lattice = (ds*tauf*Tf / (8.0*M_PI*M_PI*M_PI * norm)) * (sum_lattice + intercept_lattice);
		complex<double> result_2piDQT_05 = (ds*tauf*Tf / (8.0*M_PI*M_PI*M_PI * norm)) * (sum_2piDQT_05 + intercept_2piDQT_05);
		complex<double> result_2piDQT_10 = (ds*tauf*Tf / (8.0*M_PI*M_PI*M_PI * norm)) * (sum_2piDQT_10 + intercept_2piDQT_10);
		complex<double> result_2piDQT_15 = (ds*tauf*Tf / (8.0*M_PI*M_PI*M_PI * norm)) * (sum_2piDQT_15 + intercept_2piDQT_15);
		cout << setprecision(15) << Delta_y
				<< "   " << result_lattice.real() << "   " << result_lattice.imag()
				<< "   " << result_2piDQT_05.real() << "   " << result_2piDQT_05.imag()
				<< "   " << result_2piDQT_10.real() << "   " << result_2piDQT_10.imag()
				<< "   " << result_2piDQT_15.real() << "   " << result_2piDQT_15.imag() << endl;
	}

	return 0;
}
