/*******************************************************************
*   Propagators.h
*   KPS
*
*	Author: Kareem Omar
*	kareem.omar@uah.edu
*	https://github.com/komrad36
*
*	Last updated Feb 12, 2016
*   This application is entirely my own work.
*******************************************************************/
//
// Two generic numerical integrators with a common Propagator base class.
//
// These integrators are templated on an integer N; they operate on Arrays
// from the Eigen linear algebra library, of type Eigen::Array<double, N, 1>.
// In other words, they can solve ODEs consisting of an arbitrary number of
// elements.
//
// Initialize a propagator with init(). This function accepts a std::function
// pointing to a user-defined ODE function taking (t, y) and returning y'.
//
// Then, call step() to advance t and y by one step.
//
// The integrators available are RKDP and ABM.
//
// RKDP stands for Runge-Kutta Dormand-Prince pair and is a simple adaptive
// timestep solver. It is sufficient for most purposes and is comparable
// to MATLAB(R)'s or GNU Octave's ode45() function. It's simple and fast
// for moderate accuracies. For high accuracies, however, it's not ideal,
// particularly when ODE evaluation is expensive - which it very much is
// in KPS.
//
// Thus I also implemented ABM - an Adams-Bashforth-Moulton linear multistep solver.
// This method reuses previous step results to improve accuracy, or
// accelerate computation at equivalent accuracy, compared to RKDP.
// This is comparable to MATLAB(R)'s or GNU Octave's ode113() function.
//
// Neither of these solvers is suitable for highly stiff systems, although
// at tight tolerance and mdoerate stiffness, they may still work, just
// slowly.
//
// I have also developed a stiff solver based on the numerical differentiation
// formulas - NDF - comparable to MATLAB(R)'s or GNU Octave's ode15s() function.
// It is enormous and complicated but produces excellent results for stiff systems.
//
// HOWEVER - for only slightly stiff systems, the NDFs are both less accurate
// (due to low order) and slower, because of their complexity and the numerical
// Jacobian computation. The NDFs are only faster at tight tolerances AND stiff systems,
// such that nonstiff methods end up using way too many steps to try to deal with the stiffness.
//
// The system solved by KPS is not stiff enough to make this true; consequently,
// even at tight tolerances, it is slower than ABM and even than RKDP.
//
// NDF was a lot of work, but it's not needed in KPS, it's enormous, and it
// slows compilation times, so it has been completely removed from KPS.
//
// Contact me if interested in this solver for other purposes.
//

#ifndef PROPAGATORS_H
#define PROPAGATORS_H

#include <array>
#include <algorithm>
#include <Eigen/Eigen>
#include <functional>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <vector>


template <int N>
// base propagator class to allow branchless selection of RKDP or ABM
// 'N' is the length of the vector to be integrated
// N == 13 for KPS but this is left generic for general use
class Propagator {
private:
	typedef Eigen::Array<double, N, 1> dEvecN;
public:
	// this class contains Eigen variables, which need alignment
	// to permit use of the fastest SSE SIMD instructions
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	virtual void init(std::function<dEvecN(double, const dEvecN&)> odefun, const double t0, const dEvecN& y0, const double abs_tolerance, const double rel_tolerance, const double max_step) = 0;
	
	virtual double step(double& t_out, dEvecN& y_out) = 0;
	
	virtual ~Propagator() {}
};

template <int N>
class RKDP : public Propagator<N> {
	// --- VARIABLES ---
private:
	typedef Eigen::Array<double, N, 1> dEvecN;

	const double eps = std::numeric_limits<double>::epsilon();

	// the ODE function to call for derivative evaluation
	std::function<dEvecN(double, const dEvecN&)> ode;

	// relative error tolerance between 4th and 5th order evaluations
	double rel_tol;

	// maximum step size
	double h_max;

	// ratio of absolute to relative tolerances as threshold for new step size determination
	double tol_ratio;

	// internal time
	double t;
	
	// internal step
	double h;

	// internal vector
	dEvecN y;
	
	// minimum step size below which total error will likely begin to increase
	// excessively due to truncation error
	double h_min;

	// intermediates
	dEvecN k1, k2, k3, k4, k5, k6, k7;

	// --- /VARIABLES ---


	// --- METHODS ---
public:
	// this class contains Eigen variables, which need alignment
	// to permit use of the fastest SSE SIMD instructions
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	void init(std::function<dEvecN(double, const dEvecN&)> odefun, const double t0, const dEvecN& y0, const double abs_tolerance, const double rel_tolerance, const double max_step) {
		ode = odefun;
		t = t0;
		y = y0;
		rel_tol = rel_tolerance;
		h = h_max = max_step;
		tol_ratio = abs_tolerance / rel_tolerance;

		// initial derivative eval at t0, y0
		k1 = ode(t0, y0);

		// *** this expression reoccurs several times ***
		// We want to set min step size to approximately 16 times
		// machine epsilon. However, machine epsilon should NOT
		// really be constant; the actual epsilon, i.e. the actual
		// distance to the next representable number, depends
		// on the magnitude of the relevant value. In this case,
		// t. So instead of using raw eps, we determine it exactly
		// (for IEEE doubles) based on the number of binary digits
		// occupied by t's IEEE representation.
		h_min = t ? pow(2.0, -48.0 + floor(log2(t))) : 0.0;

		// error based on worst-performing member of first eval
		double h_inv = (k1 / y.abs().max(tol_ratio)).abs().maxCoeff() / ((4.0 / 5.0) * pow(rel_tolerance, (1.0 / 5.0)));
		
		if (h * h_inv > 1.0) h = 1.0 / h_inv;

		if (h < h_min) h = h_min;
	}

	double step(double& t_out, dEvecN& y_out) {

		// improved estimate of 16 * machine eps of variable 't'
		// see init() for full explanation
		h_min = t ? pow(2.0, -48.0 + floor(log2(t))) : 0.0;

		h = std::min(h_max, std::max(h_min, h));

		double scalar_error;
		bool step_had_failures = false;
		for (;;) {
			// compute intermediates based on Butcher tableau for RKDP
			// see https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method
			k2 = ode(t + 0.2*h, y + h*(0.2*k1));
			k3 = ode(t + 0.3*h, y + h*((3.0 / 40.0)*k1 + (9.0 / 40.0)*k2));
			k4 = ode(t + 0.8*h, y + h*((44.0 / 45.0)*k1 - (56.0 / 15.0)*k2 + (32.0 / 9.0)*k3));
			k5 = ode(t + (8.0 / 9.0)*h, y + h*((19372.0 / 6561.0)*k1 - (25360.0 / 2187.0)*k2 + (64448.0 / 6561.0)*k3 - (212.0 / 729.0)*k4));
			k6 = ode(t + h, y + h*((9017.0 / 3168.0)*k1 - (355.0 / 33.0)*k2 + (46732.0 / 5247.0)*k3 + (49.0 / 176.0)*k4 - (5103.0 / 18656.0)*k5));

			// new t
			t_out = t + h;
			
			// get step size that was *actually* stepped, i.e. due to truncation error t_out != t + h exactly
			h = t_out - t;

			// new y
			y_out = y + h * ((35.0 / 384.0)*k1 + (500.0 / 1113.0)*k3 + (125.0 / 192.0)*k4 - (2187.0 / 6784.0)*k5 + (11.0 / 84.0)*k6);

			// final derivative - for error analysis AND reusable as first intermediate of next step
			k7 = ode(t_out, y_out);

			// see Solving Ordinary Differential Equations I: Nonstiff Problems By Ernst Hairer
			dEvecN raw_error{ (71.0 / 57600.0)* k1 - (71.0 / 16695.0)* k3 + (71.0 / 1920.0)* k4 - (17253.0 / 339200.0)* k5 + (22.0 / 525.0)* k6 - (1.0 / 40.0) * k7 };

			// if worst error ratio is below rel_tol, DONE with this step!
			if ((scalar_error = h * (raw_error / y.abs().max(y_out.abs()).max(tol_ratio)).abs().maxCoeff()) <= rel_tol) break;

			// if not, shrink h and try again
			// shrink based on error as in init() if this is the first failure; otherwise, just halve it
			h = step_had_failures ? std::max(h_min, 0.5*h) : std::max(h_min, h * std::max((1.0 / 10.0), (4.0 / 5.0)*pow(rel_tol / scalar_error, (1.0 / 5.0))));
			step_had_failures = true;

		}

		// if the step succeeded immediately,
		// try to increase the step size for the next iteration
		if (!step_had_failures) {
			double factor = (5.0 / 4.0)*pow(scalar_error / rel_tol, (1.0 / 5.0));
			h = factor > (1.0 / 5.0) ? h / factor : h * 5.0;
		}

		// set internal state
		t = t_out;
		y = y_out;

		// recycle k7 for next iteration
		k1 = k7;

		// done with this step!
		return t;
	}
};

template <int N>
class ABM : public Propagator<N> {
	// --- VARIABLES ---
private:
	typedef Eigen::Array<double, N, 1> dEvecN;

	const int k_max = 12;

	const double eps = std::numeric_limits<double>::epsilon();

	const std::array<double, 13> beta_coef;

	const std::array<double, 13> powers_of_two;

	std::function<dEvecN(double, const dEvecN&)> ode;
	double rel_tol;
	double h_max;
	double tol_ratio;

	double h_prev;

	int k, k_new, k_prev;

	bool initial_ramp;

	// see Butcher p. 23
	// phi is defined as the solution Y to the IVP Y'(x)=f(x, Y(x)), Y(a)=y_0
	Eigen::Matrix<double, N, 14> phi;

	Eigen::Array<double, 12, 1> psi, alpha, beta, w, v;
	Eigen::Array<double, 13, 1> sigma, g;

	double t, h_min, h, t_prev;
	dEvecN y_last, y;

	dEvecN phi_error, pred, k1;

	dEvecN wt_inv;

	// number of steps taken at the current step size
	// (no need to change coefs if this value is greater
	// than the order)
	int num_hsteps;
		
	int num_failed_steps;

	// --- /VARIABLES ---


	// --- METHODS ---
public:
	// this class contains Eigen variables, which need alignment
	// to permit use of the fastest SSE SIMD instructions
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	// This solver uses ABM orders 1 to 12.
	// It takes guidance from the (public domain) DDEABM.F
	// by L. F. Shampine and M. K. Gordon as used in the SLATEC library,
	// a variation of which later made it into the MATLAB ODE Suite.
	// 
	// Their coefficients are actually decimal approximations with only 3 sig figs - it was the 80's, I guess
	//
	// Fractional forms are used here instead
	// with thanks to Numerical Methods for Ordinary Differential Equations, by John Butcher (2003)
	//
	// Beta Coeffs are from Theorem 410D in the Butcher text with alpha(x)=1-x
	// This is a Taylor series expansion of beta(1+z) up to 12th degree,
	ABM() : beta_coef({ (1.0 / 2.0), (1.0 / 12.0), (1.0 / 24.0), (19.0 / 720.0), (3.0 / 160.0), (863.0 / 60480.0),
	(275.0 / 24192.0), (33953.0 / 3628800.0), (8183.0 / 1036800.0), (3250433.0 / 479001600.0), (4671.0/788480.0),
	(13695779093.0/2615348736000.0), (2224234463.0/475517952000.0)}),
		powers_of_two({ 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0, 1024.0, 2048.0, 4096.0, 8192.0 }) {}

	void init(std::function<dEvecN(double, const dEvecN&)> odefun, const double t0, const dEvecN& y0, const double abs_tolerance, const double rel_tolerance, const double max_step) {
		ode = odefun;
		t = t0;
		y = y0;
		rel_tol = rel_tolerance;
		h = h_max = max_step;
		tol_ratio = abs_tolerance / rel_tolerance;

		// initial derivative eval at t0, y0
		k1 = ode(t0, y0);

		// *** this expression reoccurs several times ***
		// We want to set min step size to approximately 16 times
		// machine epsilon. However, machine epsilon should NOT
		// really be constant; the actual epsilon, i.e. the actual
		// distance to the next representable number, depends
		// on the magnitude of the relevant value. In this case,
		// t. So instead of using raw eps, we determine it exactly
		// (for IEEE doubles) based on the number of binary digits
		// occupied by t's IEEE representation.
		h_min = t ? pow(2.0, -48.0 + floor(log2(t))) : 0.0;

		// error based on worst-performing member of first eval
		double h_inv = (k1 / y.abs().max(tol_ratio)).abs().maxCoeff() / (0.25 * sqrt(rel_tolerance));

		if (h * h_inv > 1.0) h = 1.0 / h_inv;

		if (h < h_min) h = h_min;

		// start at order 1
		k = 1;

		// see Butcher p. 23
		phi.setZero();
		phi.col(0) = k1;
		sigma[0] = 1.0;
		g[0] = 1.0;
		g[1] = 0.5;

		h_prev = 0.0;
		k_prev = 0;

		// where no failures have occured yet and order (k)
		// is ramping up
		initial_ramp = true;
	}

	double step(double& t_out, dEvecN& y_out) {

		// improved estimate of 16 * machine eps of variable 't'
		// see init() for full explanation
		h_min = t ? pow(2.0, -48.0 + floor(log2(t))) : 0.0;

		h = std::min(h_max, std::max(h_min, h));

		num_failed_steps = 0;
		wt_inv = y.abs().max(tol_ratio).cwiseInverse();
		double k_minus_2_error, k_minus_1_error, k_error, k_plus_1_error;
		for (;;) {
			if (h != h_prev) num_hsteps = 0;
			if (num_hsteps <= k_prev) ++num_hsteps;

			// update coefs
			if (k >= num_hsteps) {
				beta[num_hsteps - 1] = 1.0;
				alpha[num_hsteps - 1] = 1.0 / num_hsteps;
				double h_times_num_hsteps = h * num_hsteps;
				sigma[num_hsteps] = 1.0;
				for (int i = num_hsteps + 1; i <= k; ++i) {
					double old_psi_i_minus_2 = psi[i - 2];
					psi[i - 2] = h_times_num_hsteps;
					h_times_num_hsteps = old_psi_i_minus_2 + h;

					beta[i - 1] = beta[i - 2] * psi[i - 2] / old_psi_i_minus_2;
					alpha[i - 1] = h / h_times_num_hsteps;
					sigma[i] = i * alpha[i - 1] * sigma[i - 1];
				}
				psi[k - 1] = h_times_num_hsteps;

				// new g coefs from v and w
				if (num_hsteps == 1) {
					for (int i = 1; i <= k; ++i) {
						v[i - 1] = 1.0 / static_cast<double>(i * (i + 1));
					}
					w = v;
				}
				else {
					// diag v only changes when order increases
					if (k > k_prev) {
						v[k - 1] = 1.0 / static_cast<double>(k * (k + 1));
						for (int j = 1; j <= num_hsteps - 2; ++j) {
							v[k - j - 1] -= alpha[j] * v[k - j];
						}
					}
					for (int iq = 1; iq <= k + 1 - num_hsteps; ++iq) {
						w[iq - 1] = v[iq - 1] -= alpha[num_hsteps - 1] * v[iq];
					}
					g[num_hsteps] = w[0];
				}

				// new g from working vector
				for (int i = num_hsteps + 2; i <= k + 1; ++i) {
					for (int iq = 1; iq <= k + 2 - i; ++iq) {
						w[iq - 1] -= alpha[i - 2] * w[iq];
					}
					g[i - 1] = w[0];
				}
			}

			// updated phi from betas
			for (int i = num_hsteps; i < k; ++i) {
				phi.col(i) *= beta[i];
			}

			// --- PE STEP ---
			// phi solution PE
			phi.col(k + 1) = phi.col(k);
			phi.col(k).setZero();
			pred.setZero();

			// accumulate from phi
			for (int i = k; i >= 1; --i) {
				pred += g[i - 1] * phi.col(i - 1).array();
				phi.col(i - 1) += phi.col(i);
			}

			pred = y + h * pred;

			t_prev = t;
			t = t_prev + h;

			k1 = ode(t, pred);

			// errors from k-2 to k
			phi_error = k1 - phi.col(0).array();
			double wt_phi_max = (wt_inv*phi_error).abs().maxCoeff();
			double scalar_error = h * (g[k - 1] - g[k]) * wt_phi_max;
			k_error = h * sigma[k] * beta_coef[k - 1] * wt_phi_max;

			k_minus_1_error = k >= 2 ? (h * sigma[k - 1] * beta_coef[k - 2] * ((phi.col(k - 1).array() + phi_error) * wt_inv).abs().maxCoeff()) : 0.0;
			k_minus_2_error = k >= 3 ? (h * sigma[k - 2] * beta_coef[k - 3] * ((phi.col(k - 2).array() + phi_error) * wt_inv).abs().maxCoeff()) : 0.0;

			// need to lower error?
			k_new = k - (((k == 2) && (k_minus_1_error <= 0.5*k_error)) || ((k > 2) && std::max(k_minus_1_error, k_minus_2_error) <= k_error));

			// if error at order k is below rel_tol, DONE with this step!
			if (scalar_error <= rel_tol) break;

			// otherwise, have to roll back t, phi, psi, and shrink step
			// and, if enough failures, reset the order
			initial_ramp = false;
			t = t_prev;
			for (int i = 1; i <= k; ++i) {
				phi.col(i - 1) = (phi.col(i - 1) - phi.col(i)) / beta[i - 1];
			}
			for (int i = 2; i <= k; ++i) {
				psi[i - 2] = psi[i - 1] - h;
			}

			double shrink_factor = 0.5;
			if (++num_failed_steps == 3) {
				k_new = 1;
			}
			else if (num_failed_steps > 3) {
				shrink_factor = std::min(0.5, sqrt(0.5*rel_tol / k_error));
			}
			h = std::max(shrink_factor * h, h_min);
			k = k_new;
		}

		k_prev = k;
		h_prev = h;
		y_last = y;

		// --- CE STEP ---
		y = pred + h * g[k] * phi_error;

		// improved ODE eval
		k1 = ode(t, y);

		// updated phi sol diffs
		phi.col(k) = k1 - phi.col(0).array();
		phi.col(k + 1) = phi.col(k) - phi.col(k + 1);
		for (int i = 1; i <= k; ++i) {
			phi.col(i - 1) += phi.col(k);
		}

		// if we lost an order or are at max order, don't raise order
		// and exit initial ramp stage
		if (k_new == k - 1 || k == k_max) initial_ramp = false;

		// determine new order
		if (initial_ramp) {
			++k;
		}
		else if (k_new == k - 1) {
			--k;
			k_error = k_minus_1_error;
		}
		else if (k + 1 <= num_hsteps) {
			k_plus_1_error = h * beta_coef[k] * (phi.col(k + 1).array() * wt_inv).abs().maxCoeff();
			if (k == 1) {
				if (k_plus_1_error < 0.5 * k_error) {
					++k;
					k_error = k_plus_1_error;
				}
			}
			else {
				if (k_minus_1_error <= std::min(k_error, k_plus_1_error)) {
					--k;
					k_error = k_minus_1_error;
				}
				else if (k < k_max && k_plus_1_error < k_error) {
					++k;
					k_error = k_plus_1_error;
				}
			}
		}

		// determine new step
		// double h if in initial stage or low error
		if (initial_ramp || (k_error * powers_of_two[k] <= 0.5 * rel_tol)) {
			h *= 2.0;
		}
		// shrink h if high error
		else if (k_error > 0.5*rel_tol) {
			h *= std::max(0.5, std::min((9.0 / 10.0), pow(0.5 * rel_tol / k_error, 1.0 / (k + 1))));
		}

		// report new t and y
		t_out = t;
		y_out = y;

		// done with this step!
		return t;
	}
};

#endif
