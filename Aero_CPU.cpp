/*******************************************************************
*   Aero_CPU.cpp
*   KPS
*
*	Author: Kareem Omar
*	kareem.omar@uah.edu
*	https://github.com/komrad36
*
*	Last updated Feb 27, 2016
*   This application is entirely my own work.
*******************************************************************/
//
// Aero uses CUDA or CPU and operates on a satellite or body defined as a series of polygons
// in 3-D space. When supplied with air density and velocity (in the Body frame),
// it approximates the drag force and torque on the body by simulating
// collisions and accumulating impulses and angular impulses per unit time.
//
// Note that net force and torque are returned in the Body frame.
//
// This module is intended for use by a 6 DoF orbital/attitude propagator which
// calls the aer() function from its integrators to obtain forces and/or torques,
// such as KPS.
//

#include "Aero.h"

#define PAD				(0.000001)

// placeholder x-value to signal no collision occurred
#define NONE			(-9999999.0)

// internal collision routine
void Aero_CPU::collide(vec3& f, vec3& t, const double rho, const double v_mag2) {
	vec3 best_N, r;

	vec3 sum_F_compensation{ 0.0, 0.0, 0.0 };
	vec3 force, interm_F, v_norm;

	vec3 sum_T_compensation{ 0.0, 0.0, 0.0 };
	vec3 torque, interm_T;

	// for each particle within the bounding box of the satellite...
	for (int k = min_k; k <= max_k; ++k) {
		double y = k * pitch;
		for (int m = min_m; m <= max_m; ++m) {
			double z = m * pitch;
			double best_x = NONE;
			// ...check for collisions against each poly:
			for (int p = 0; p < num_poly; ++p) {
				// bail early if outside bounding box of this particular poly...
				if (y < min_y[p] || y > max_y[p] || z < min_z[p] || z > max_z[p]) continue;

				// ...otherwise, perform point-in-polygon anaylsis.
				// Looks scary but isn't too bad. Pretend coordinate system
				// is such that test point is at origin. Send the point to the right (+x)
				// and flip a flag each time it crosses a segment of the polygon. Flag will
				// end up flipped if inside polygon because it will cross an odd number of
				// segments. This works even for concave polygons.
				// For each segment of the polygon:
				//		if both ends are above or below x axis (i.e. share same sign), not a crossing
				//			otherwise, if both ends are left of the y axis, not a crossing
				//				otherwise, we do the slow bit (but usually don't have to due to the above):
				//				see if the segment intersects the x axis right of 0. if it does, crossing!
				// See KPS research paper for more.
				vec3* P = P_rot + p*NUM_VTX;
				int j = NUM_VTX - 1;
				int odd_nodes = 0;
				for (int i = 0; i < NUM_VTX; ++i) {
					odd_nodes ^= (((((P[i].z < z && P[j].z >= z) || (P[j].z < z && P[i].z >= z))
						&& (P[i].y <= y || P[j].y <= y)))
						&& ((P[i].y + (z - P[i].z) / (P[j].z - P[i].z)*(P[j].y - P[i].y) < y)));
					j = i;
				}

				// if inside polygon, compute collision location
				if (odd_nodes && N[p].x) {
					double x = (precomp[p] - (N[p].y*y + N[p].z*z)) / N[p].x;

					// and if it's the best (first) collision so far, update best_x
					if (x > best_x) {
						best_x = x;
						best_N = N[p];
					}
				}
			}

			// if a collision occurred
			if (best_x > NONE + PAD) {
				// Kahan sum of force to minimize numerical error
				// see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				//
				// see Equation 62 in KPS Research paper
				force = -f_scalar*rho*v_mag2*best_N.x*best_N - sum_F_compensation;
				interm_F = f + force;
				sum_F_compensation = (interm_F - f) - force;
				f = interm_F;

				// Kahan sum of torque to minimize numerical error
				// see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
				//
				// see Equation 63 in KPS Research paper
				torque = glm::cross(vec3{ best_x, y, z } - CM_R, force) - sum_T_compensation;
				interm_T = t + torque;
				sum_T_compensation = (interm_T - t) - torque;
				t = interm_T;
			}
		}
	}
}

// collide simulated particles to approximate the resulting forces and torques,
// 'f' and 't'
void Aero_CPU::aer(vec3& f, vec3& t, const double rho, const vec3& v) {

	// zero out 'f' and 't'
	f = t = vec3();

	double v_mag2 = glm::length2(v);
	double v_mag = sqrt(v_mag2);

	// --- ROTATOR SETUP ---
	// rotates entire satellite such that velocity is in +x in new frame, i.e. relative wind
	// direction is -x
	//
	// Uses Rodrigues' rotation formula for speed
	// the ternary handles cases where velocity is entirely in +x or entirely in -x already
	//
	// does NOT handle zero velocity, but if that's the case in an *orbital simulation*
	// you have bigger problems to worry about
	double cos_theta = v.x / v_mag;
	double sin_theta = sin(acos(cos_theta));
	vec3 k = sin_theta ? vec3{0.0, v.z, -v.y} / (v_mag * sin_theta) : vec3();
	vec3 k_times_1_minus_cos_theta = k*(1 - cos_theta);
	// --- /ROTATOR SETUP ---

	// Use of ROTATOR (Rodrigues' formula) to rotate satellite CM into new frame
	CM_R = cos_theta*CM + sin_theta*glm::cross(k, CM) + k_times_1_minus_cos_theta*glm::dot(k, CM);

	for (int i = 0; i < total_pts; ++i) {
		// Use of ROTATOR (Rodrigues' formula) to rotate polygons into new frame
		P_rot[i] = cos_theta*P_s[i] + sin_theta*glm::cross(k, P_s[i]) + k_times_1_minus_cos_theta*glm::dot(k, P_s[i]);
	}

	for (int i = 0; i < num_poly; ++i) {
		precompute(i);
	}

	// get whole satellite's bounding box
	double total_min_y = *std::min_element(min_y, min_y + num_poly) + PAD;
	double total_max_y = *std::max_element(max_y, max_y + num_poly) - PAD;
	double total_min_z = *std::min_element(min_z, min_z + num_poly) + PAD;
	double total_max_z = *std::max_element(max_z, max_z + num_poly) - PAD;

	// determine number of simulated particles required in each direction based on linear pitch
	min_k = static_cast<int>(total_min_y / pitch);
	max_k = static_cast<int>(total_max_y / pitch);
	min_m = static_cast<int>(total_min_z / pitch);
	max_m = static_cast<int>(total_max_z / pitch);

	collide(f, t, rho, v_mag2);

	// Use of ROTATOR (Rodrigues' formula) to rotate summed force and torque BACK TO BODY FRAME
	// (note the negative sign)
	f = cos_theta*f - sin_theta*glm::cross(k, f) + k_times_1_minus_cos_theta*glm::dot(k, f);
	t = cos_theta*t - sin_theta*glm::cross(k, t) + k_times_1_minus_cos_theta*glm::dot(k, t);
}