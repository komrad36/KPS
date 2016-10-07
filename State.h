/*******************************************************************
*   State.h
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
// Dual representation of a satellite state as:
//		a)	a Plain Old Data vector of 13 doubles,
//			i.e. an Eigen Array or Matrix type with dimensions
//			13 by 1, as needed by the numerical integrators, which
//			do not care about the meaning of the different
//			parts of the vector and integrate it all in one go
//		b)	an organized Plain Old Data struct consisting of
//			the satellite's position, velocity, orientation,
//			and angular velocity, represented by glm vectors
//			and a glm quaternion (for the orientation)
//		
// The ability to reinterpret the same struct without copying
// is crucial to the speed of the application. Parts of the
// 13-vector, such as the velocity, say, can be rapidly treated as
// individual 3-vectors for manipulation by the ODE, printing, etc.,
// yet when Eigen needs to consider the whole vector for integration,
// it can do so with no overhead.
//

#pragma once

#include <algorithm>

#include "glm_util.h"
#include "Eigen_util.h"

// Eigen requires 16-byte alignment for SSE. Because this
// struct will be reinterpreted as an Eigen Array, it, too,
// must align itself along 16-byte boundaries.
//
// NOTE: due to this pack request the actual State struct will
// end up 14 doubles long, not 13. This is okay. The padding is at
// the end and does not interfere with anything in conversion.
struct alignas(16) State {
	// Position in ECI frame [m]
	vec3 r;

	// Velocity in ECI frame [m/s]
	vec3 v;

	// Orientation
	// NOTE: this quaternion turns the ECI frame axes into the Body frame axes,
	// i.e. rotates a VECTOR from the Body frame into the ECI frame.
	// See the KPS research paper for more.
	quat q;

	// Angular velocity of the satellite about the body axes [rad/s]
	// i.e. angular velocity of Body frame relative to ECI frame,
	// coordinatized in the body frame
	vec3 w;

	State() {}
	State(const vec3& state_r, const vec3& state_v, const quat& state_q, const vec3& state_w) : r(state_r), v(state_v), q(state_q), w(state_w) {}
};

// Interpret a const Eigen 13-vector as a const State struct
inline const State& e_to_s(const dEvec13& e_state) {
	return reinterpret_cast<const State&>(e_state);
}

// Interpret an Eigen 13-vector as a State struct
inline State& e_to_s(dEvec13& e_state) {
	return reinterpret_cast<State&>(e_state);
}

// Interpret a const State struct as a const Eigen 13-vector
inline const dEvec13& s_to_e(const State& s_state) {
	return reinterpret_cast<const dEvec13&>(s_state);
}

// Interpret a State struct as an Eigen 13-vector
inline dEvec13& s_to_e(State& s_state) {
	return reinterpret_cast<dEvec13&>(s_state);
}

