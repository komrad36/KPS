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

// collide simulated particles to approximate the resulting forces and torques,
// 'f' and 't'
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

	ClipperLib::Path clip, subj;
	clip.reserve(NUM_VTX);
	subj.reserve(NUM_VTX);
	ClipperLib::Paths clip_list;
	clip_list.reserve(num_poly);
	// for each polygon
	for (int i = 0; i < num_poly; ++i) {
		// check other polys and build a list of those in front of ith polygon

		// subject polygon
		vec3* S = P_rot + i*NUM_VTX;
		subj.clear();
		clip_list.clear();
		for (int vtx = 0; vtx < NUM_VTX; ++vtx) {
			subj.push_back(ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*S[vtx].y), static_cast<ClipperLib::cInt>(CLIP_MUL*S[vtx].z)));
		}

		for (int j = (i == 0); j < num_poly; j += j==i-1 ? 2 : 1) {

			// candidate polygon
			vec3* C = P_rot + j*NUM_VTX;

			if (j < i) {
				if (!lower_idx_is_above[j*num_poly + i]) {
					clip_list.push_back({
						ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*C[0].y), static_cast<ClipperLib::cInt>(CLIP_MUL*C[0].z)),
						ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*C[1].y), static_cast<ClipperLib::cInt>(CLIP_MUL*C[1].z)),
						ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*C[2].y), static_cast<ClipperLib::cInt>(CLIP_MUL*C[2].z)),
						ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*C[3].y), static_cast<ClipperLib::cInt>(CLIP_MUL*C[3].z))
					});
				}
				continue;
			}

			lower_idx_is_above[i*num_poly + j] = true;

			// if no overlap from this candidate, do not add to clip list
			if (min_y[i] >= max_y[j] || max_y[i] <= min_y[j] || min_z[i] >= max_z[j] || max_z[i] <= min_z[j]) continue;

			// if there is overlap, check vertices of intersection to see which polygon is above.
			// if candidate is below, do not add to clip list
			// if candidate is above, add to clip list
			clip = {
				ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*C[0].y), static_cast<ClipperLib::cInt>(CLIP_MUL*C[0].z)),
				ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*C[1].y), static_cast<ClipperLib::cInt>(CLIP_MUL*C[1].z)),
				ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*C[2].y), static_cast<ClipperLib::cInt>(CLIP_MUL*C[2].z)),
				ClipperLib::IntPoint(static_cast<ClipperLib::cInt>(CLIP_MUL*C[3].y), static_cast<ClipperLib::cInt>(CLIP_MUL*C[3].z))
			};
			clipper.Clear();
			ClipperLib::Paths sol;
			clipper.AddPath(subj, ClipperLib::ptSubject, true);
			clipper.AddPath(clip, ClipperLib::ptClip, true);
			clipper.Execute(ClipperLib::ctIntersection, sol);
			if (sol.empty()) continue;
			lower_idx_is_above[i*num_poly + j] = false;
			for (const ClipperLib::IntPoint& pt : sol[0]) {
				double y = INV_CLIP_MUL*static_cast<double>(pt.X);
				double z = INV_CLIP_MUL*static_cast<double>(pt.Y);
				double subj_x = (precomp[i] - (N[i].y*y + N[i].z*z)) / N[i].x;
				double clip_x = (precomp[j] - (N[j].y*y + N[j].z*z)) / N[j].x;
				if (subj_x > clip_x) break;
				if (subj_x < clip_x) {
					lower_idx_is_above[i*num_poly + j] = true;
					clip_list.push_back(clip);
					break;
				}
			}
		}

		// now that clip list is complete, get unoccluded portion
		clipper.Clear();
		ClipperLib::Paths sol;
		clipper.AddPath(subj, ClipperLib::ptSubject, true);
		clipper.AddPaths(clip_list, ClipperLib::ptClip, true);
		clipper.Execute(ClipperLib::ctDifference, sol);
		for (const ClipperLib::Path& path : sol) {
			ClipperLib::cInt area = 0, cm_y = 0, cm_z = 0;
			size_t highest = path.size() - 1;
			size_t m = highest;
			for (size_t n = 0; n <= highest; ++n) {
				ClipperLib::cInt factor = path[n].X*path[m].Y - path[m].X*path[n].Y;
				area += factor;
				cm_y += (path[n].X + path[m].X)*factor;
				cm_z += (path[n].Y + path[m].Y)*factor;
				m = n;
			}
			double area_d = INV_CLIP_MUL*INV_CLIP_MUL*0.5*area;
			double cm_d_y = INV_CLIP_MUL*INV_CLIP_MUL*INV_CLIP_MUL*static_cast<double>(cm_y) / (6.0 * area_d);
			double cm_d_z = INV_CLIP_MUL*INV_CLIP_MUL*INV_CLIP_MUL*static_cast<double>(cm_z) / (6.0 * area_d);

			vec3 force = 2.0*fabs(area_d)*rho*N[i] * (-v_mag2*N[i].x*fabs(N[i].x));
			f += force;
			if (N[i].x) t += glm::cross(vec3{ (precomp[i] - (N[i].y*cm_d_y + N[i].z*cm_d_z)) / N[i].x, cm_d_y, cm_d_z } -CM_R, force);
		}

	}

	// --- /COMPUTE ---


	// Use of ROTATOR (Rodrigues' formula) to rotate summed force and torque BACK TO BODY FRAME
	// (note the negative sign)
	f = cos_theta*f - sin_theta*glm::cross(k, f) + k_times_1_minus_cos_theta*glm::dot(k, f);
	t = cos_theta*t - sin_theta*glm::cross(k, t) + k_times_1_minus_cos_theta*glm::dot(k, t);

}
