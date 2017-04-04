/*******************************************************************
*   Aero_Analytical.cpp
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
// Aero uses collision-based or analytical modes and operates on a satellite
// or body defined as a series of polygons
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

// analyze visible satellite geometry to determine the resulting forces and torques,
// 'f' and 't', assuming specular reflection
void Aero_Analytical::aer(vec3& f, vec3& t, const double rho, const vec3& v) {

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
	vec3 k = sin_theta ? vec3{ 0.0, v.z, -v.y } / (v_mag * sin_theta) : vec3();
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

	// --- COMPUTE ---

	// zero out 'f' and 't'
	f = t = vec3();

	ClipperLib::Paths sol;
	// consider each polygon in turn, referred to as the 'subject' polygon.
	for (int i = 0; i < num_poly; ++i) {

		// subject polygon ptr
		vec3* S = P_rot + i*NUM_VTX;

		// convert to integer coords as required by Clipper library to prevent numerical error
		// currently hardcoded to 4 vertices for speed; modify if you change NUM_VTX,
		// or just substitute a for loop to handle all cases automatically (slower)
		// or a template solution (more complicated)
		//
		// CLIP_MUL is chosen to minimize numerical error and remain within the bounds of 64-bit signed integers
		// for typical polygon coordinates
		ClipperLib::Paths subjs = { {
				ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*S[0].y), static_cast<ClipperLib::cInt>(CLIP_MUL*S[0].z)),
				ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*S[1].y), static_cast<ClipperLib::cInt>(CLIP_MUL*S[1].z)),
				ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*S[2].y), static_cast<ClipperLib::cInt>(CLIP_MUL*S[2].z)),
				ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*S[3].y), static_cast<ClipperLib::cInt>(CLIP_MUL*S[3].z))
			} };

		// now consider each of the remaining polygons in turn, referred to as the 'candidate' polygon
		for (int j = (i == 0); j < num_poly; j += 1 + (j==i-1)) {

			// candidate polygon ptr
			vec3* C = P_rot + j*NUM_VTX;

			// if j < i we already have information on which polygon is in front
			// and, thus, whether the subject polygon can be occluded by this candidate
			if (j < i) {
				if (!lower_idx_is_above[j*num_poly + i]) {
					clipper.Clear();
					clipper.AddPaths(subjs, ClipperLib::ptSubject, true);
					clipper.AddPath({
						ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*C[0].y), static_cast<ClipperLib::cInt>(CLIP_MUL*C[0].z)),
						ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*C[1].y), static_cast<ClipperLib::cInt>(CLIP_MUL*C[1].z)),
						ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*C[2].y), static_cast<ClipperLib::cInt>(CLIP_MUL*C[2].z)),
						ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*C[3].y), static_cast<ClipperLib::cInt>(CLIP_MUL*C[3].z))
					}, ClipperLib::ptClip, true);
					clipper.Execute(ClipperLib::ctDifference, subjs);
				}
			}
			else {
				lower_idx_is_above[i*num_poly + j] = true;

				// if the bounding boxes of subject and candidate don't overlap, there's no occlusion,
				// so we can bail out early, which improves performance despite the cost of evaluating the bounding extrema
				if (min_y[i] >= max_y[j] || max_y[i] <= min_y[j] || min_z[i] >= max_z[j] || max_z[i] <= min_z[j]) continue;

				// otherwise, use Clipper to determine the intersection of the two polygons,
				// and check the intersection's vertices to see which polygon is in front
				clipper.Clear();
				clipper.AddPaths(subjs, ClipperLib::ptSubject, true);
				clipper.AddPath({
					ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*C[0].y), static_cast<ClipperLib::cInt>(CLIP_MUL*C[0].z)),
					ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*C[1].y), static_cast<ClipperLib::cInt>(CLIP_MUL*C[1].z)),
					ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*C[2].y), static_cast<ClipperLib::cInt>(CLIP_MUL*C[2].z)),
					ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*C[3].y), static_cast<ClipperLib::cInt>(CLIP_MUL*C[3].z))
				},ClipperLib::ptClip, true);
				clipper.Execute(ClipperLib::ctIntersection, sol);
				if (sol.empty()) continue;
				for (auto&& pt : sol[0]) {
					double y = INV_CLIP_MUL*static_cast<double>(pt.X);
					double z = INV_CLIP_MUL*static_cast<double>(pt.Y);
					double subj_x = (precomp[i] - (N[i].y*y + N[i].z*z)) / N[i].x;
					double clip_x = (precomp[j] - (N[j].y*y + N[j].z*z)) / N[j].x;

					// if the subject is in front, the candidate does not occlude it,
					// but the candidate does occlude the subject, so store that info
					// to save time later
					//
					// numerical error is once again present in this floating point
					// computation so CLIP_PAD prevents false triggers when the
					// points are actually equal (so we should do nothing and wait
					// for another vertex to be tested that is not shared between
					// subject and candidate), but end up very slightly apart.
					if (subj_x > clip_x + CLIP_PAD) {
						lower_idx_is_above[i*num_poly + j] = false;
						break;
					}
					// if the candidate is in front, clip the subject with the
					// candidate to get the remaining unoccluded polygon 
					if (subj_x + CLIP_PAD < clip_x) {
						clipper.Execute(ClipperLib::ctDifference, subjs);
						break;
					}
				}
			}
		}

		// subjs now has one or more polygons representing the portion of the original
		// subject polygon that remain unoccluded by any other polygons.
		// 
		// thus, subjs contains panels exposed to aerodynamic force.
		// for each one, compute the area and center of mass,
		// then use that information to compute total force and torque on that panel
		// and accumulate the totals in 'f' and 't'
		for (auto&& path : subjs) {
			double area = 0.0, cm_y = 0.0, cm_z = 0.0;
			size_t highest = path.size() - 1;
			size_t m = highest;
			for (size_t n = 0; n <= highest; ++n) {
				double x_n = INV_CLIP_MUL*static_cast<double>(path[n].X);
				double x_m = INV_CLIP_MUL*static_cast<double>(path[m].X);
				double y_n = INV_CLIP_MUL*static_cast<double>(path[n].Y);
				double y_m = INV_CLIP_MUL*static_cast<double>(path[m].Y);

				double factor = x_n*y_m - x_m*y_n;
				area += factor;
				cm_y += (x_n + x_m)*factor;
				cm_z += (y_n + y_m)*factor;

				m = n;
			}
			area *= 0.5;
			double inv_6_area = area ? 1.0 / (6.0 * area) : 0.0;
			cm_y *= inv_6_area;
			cm_z *= inv_6_area;

			vec3 force = -2.0*rho*v_mag2*fabs(area)*N[i].x*N[i];
			f += force;
			if (N[i].x) t += glm::cross(vec3{ (precomp[i] - (N[i].y*cm_y + N[i].z*cm_z)) / N[i].x, cm_y, cm_z } -CM_R, force);
		}
	}

	// --- /COMPUTE ---

	// Use of ROTATOR (Rodrigues' formula) to rotate summed force and torque BACK TO BODY FRAME
	// (note the negative sign)
	f = cos_theta*f - sin_theta*glm::cross(k, f) + k_times_1_minus_cos_theta*glm::dot(k, f);
	t = cos_theta*t - sin_theta*glm::cross(k, t) + k_times_1_minus_cos_theta*glm::dot(k, t);
}
