/*******************************************************************
*   Aero.h
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
// Aero can also analytically compute these results by obtaining non-occluded frontal areas
// using the Clipper polygon clipping library.
//
// Note that net force and torque are returned in the Body frame.
//
// This module is intended for use by a 6 DoF orbital/attitude propagator which
// calls the aer() function from its integrators to obtain forces and/or torques,
// such as KPS.
//

#pragma once

#define CLIP_MUL 1e6
#define INV_CLIP_MUL 1e-6

#include "glm_util.h"

#include "clipper.hpp"

#include <iostream>
#include <thread>
#include <algorithm>
#include <iostream>
#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

// compute the next highest power of 2 above 'x'
// returns 'x' if it's already a power of 2
//
// useful for CUDA when it's necessary to dispatch
// a power-of-two number of threads (for reductions)
inline int32_t nextHigherPow2(int32_t x) {
	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return ++x;
}

// pure virtual base class allowing branchless selection of CPU vs CUDA vs Analytical
class Aero {
public:
	virtual void aer(vec3& f, vec3& t, const double rho, const vec3& v) = 0;
	virtual ~Aero() {}
};

class Aero_CUDA : public Aero {
	// --- VARIABLES ---
private:

	// device ptr to rotated polygons
	vec3* d_P_rot;

	// device ptr to forces accumulated by a CUDA block
	vec3* d_block_sums_f;

	// host ptr to forces accumulated by a CUDA block
	vec3* block_sums_f;

	// device ptr to torques accumulated by a CUDA block
	vec3* d_block_sums_t;

	// host ptr to torques accumulated by a CUDA block
	vec3* block_sums_t;

	// device ptr to polygon normal vectors
	vec3* d_N;

	// device ptr to precomputation results
	double* d_precomp;

	// device ptrs to polygon bounds
	double* d_min_y;
	double* d_max_y;
	double* d_min_z;
	double* d_max_z;

	// device ptr to accumulated total bounds
	double* d_totals;

	// host array of accumulated total bounds
	double totals[4];

	// device ptr to rotated satellite center of mass
	vec3* d_CM_R;

	// host ptr to Body frame polygons
	vec3* const P_s;

	// host ptr to rotated polygons
	vec3* const P_rot;

	// CUDA Device ID
	int cuda_device;

	// linear pitch between collisions
	const double pitch;

	// scalar constant (2A = 2*pitch^2) that multiplies force
	// see KPS research paper
	const double f_scalar;

	const int num_poly;

	const int total_pts;

	// host storage of satellite center of mass
	vec3 CM;

	// --- /VARIABLES ---


	// --- METHODS ---
public:
	Aero_CUDA(const double linear_pitch, const int num_polygons, const vec3* const poly, const vec3 sat_CM);

	~Aero_CUDA();

	void aer(vec3& f, vec3& t, const double rho, const vec3& v);

	bool init(const int cuda_device);

	// --- /METHODS ---
};

// functionoid to compare two vectors by their y-component
struct compare_y {
	bool operator()(vec3& lhs, vec3& rhs) {
		return lhs.y < rhs.y;
	}
};

// functionoid to compare two vectors by their z-component
struct compare_z{
	bool operator()(vec3& lhs, vec3& rhs) {
		return lhs.z < rhs.z;
	}
};

class Aero_CPU : public Aero {
	// --- VARIABLES ---
private:

	// loop step bounds determined by polygon bounds
	int min_k;
	int max_k;
	int min_m;
	int max_m;

	// rotated satellite center of mass
	vec3 CM_R;

	// linear pitch between collisions
	const double pitch;

	const int num_poly;

	// scalar constant (2A = 2*pitch^2) that multiplies force
	// see KPS research paper
	const double f_scalar;

	const int total_pts;

	// satellite center of mass
	const vec3 CM;

	// ptr to Body frame polygons
	vec3* const P_s;

	// ptr to polygon normal vectors
	vec3* const N;

	// ptr to precomputation results
	double* const precomp;

	// ptrs to polygon bounds
	double* const min_y;
	double* const max_y;
	double* const min_z;
	double* const max_z;

	// ptr to rotated polygons
	vec3* const P_rot;

	// --- /VARIABLES ---


	// --- METHODS ---
private:
	// precompute some info for speed,
	// including panel normals and some of collision location math
	inline void precompute(const int i) {
		vec3* P = P_rot + i*NUM_VTX;
		N[i] = glm::normalize(glm::cross(P[1] - P[0], P[1] - P[2]));
		precomp[i] = glm::dot(N[i], P[0]);

		// find mins and maxes of y and z of each panel
		std::pair<vec3*, vec3*> y_pair = std::minmax_element(P, P + NUM_VTX, compare_y());
		std::pair<vec3*, vec3*> z_pair = std::minmax_element(P, P + NUM_VTX, compare_z());
		min_y[i] = y_pair.first->y;
		max_y[i] = y_pair.second->y;
		min_z[i] = z_pair.first->z;
		max_z[i] = z_pair.second->z;
	}

	void collide(vec3& f, vec3& t, const double rho, const double v_mag2);

public:

	Aero_CPU(const double linear_pitch, const int num_polygons, const vec3* const poly, const vec3 sat_CM) :
		pitch(linear_pitch),
		num_poly(num_polygons),
		f_scalar(2.0*linear_pitch*linear_pitch),
		total_pts(num_polygons * NUM_VTX),
		CM(sat_CM),
		P_s(new vec3[num_polygons * NUM_VTX]),
		N(new vec3[num_polygons]),
		precomp(new double[num_polygons]),
		min_y(new double[num_polygons]),
		max_y(new double[num_polygons]),
		min_z(new double[num_polygons]),
		max_z(new double[num_polygons]),
		P_rot(new vec3[num_polygons * NUM_VTX]) {

		// copy polygon data into internal storage
		for (int i = 0; i < total_pts; ++i) {
			P_s[i] = poly[i];
		}
	}

	~Aero_CPU() {
		delete[] P_s;
		delete[] N;
		delete[] precomp;
		delete[] min_y;
		delete[] max_y;
		delete[] min_z;
		delete[] max_z;
		delete[] P_rot;
	}

	void aer(vec3& f, vec3& t, const double rho, const vec3& v);

	// --- /METHODS ---

};


class Aero_Analytical : public Aero {
	// --- VARIABLES ---
private:

	// rotated satellite center of mass
	vec3 CM_R;

	const int num_poly;

	const int total_pts;

	// satellite center of mass
	const vec3 CM;

	// ptr to Body frame polygons
	vec3* const P_s;

	// ptr to polygon normal vectors
	vec3* const N;

	// ptr to precomputation results
	double* const precomp;

	// ptrs to polygon bounds
	double* const min_y;
	double* const max_y;
	double* const min_z;
	double* const max_z;

	// ptr to rotated polygons
	vec3* const P_rot;

	bool* const lower_idx_is_above;

	ClipperLib::Clipper clipper;

	// --- /VARIABLES ---


	// --- METHODS ---
private:
	// precompute some info for speed,
	// including panel normals and extents
	inline void precompute(const int i) {
		vec3* P = P_rot + i*NUM_VTX;
		N[i] = glm::normalize(glm::cross(P[1] - P[0], P[1] - P[2]));
		precomp[i] = glm::dot(N[i], P[0]);

		// find mins and maxes of y and z of each panel
		std::pair<vec3*, vec3*> y_pair = std::minmax_element(P, P + NUM_VTX, compare_y());
		std::pair<vec3*, vec3*> z_pair = std::minmax_element(P, P + NUM_VTX, compare_z());
		min_y[i] = y_pair.first->y;
		max_y[i] = y_pair.second->y;
		min_z[i] = z_pair.first->z;
		max_z[i] = z_pair.second->z;
	}

public:

	Aero_Analytical(const double linear_pitch, const int num_polygons, const vec3* const poly, const vec3 sat_CM) :
		num_poly(num_polygons),
		total_pts(num_polygons * NUM_VTX),
		CM(sat_CM),
		P_s(new vec3[num_polygons * NUM_VTX]),
		N(new vec3[num_polygons]),
		precomp(new double[num_polygons]),
		min_y(new double[num_polygons]),
		max_y(new double[num_polygons]),
		min_z(new double[num_polygons]),
		max_z(new double[num_polygons]),
		P_rot(new vec3[num_polygons * NUM_VTX]),
		lower_idx_is_above(new bool[num_polygons * num_polygons]) {

		// silence unused parameter
		static_cast<void>(linear_pitch);
	
		// copy polygon data into internal storage
		for (int i = 0; i < total_pts; ++i) {
			P_s[i] = poly[i];
		}
	}

	~Aero_Analytical() {
		delete[] P_s;
		delete[] N;
		delete[] precomp;
		delete[] min_y;
		delete[] max_y;
		delete[] min_z;
		delete[] max_z;
		delete[] P_rot;
		delete[] lower_idx_is_above;
	}

	void aer(vec3& f, vec3& t, const double rho, const vec3& v);

	// --- /METHODS ---

};
