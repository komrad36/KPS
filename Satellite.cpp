/*******************************************************************
*   Satellite.cpp
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
// Satellite model, including satellite orbital and attitude ODE.
//
 
#include "Satellite.h"

// check if the satellite is still in orbit
bool Satellite::isInValidState() const {
	const double deorbit_r = earth.DEORBIT_ALT + earth.R;
	return glm::length2(state.r) > deorbit_r*deorbit_r;
}

// the derivative of an Eigen 13-vector representing the full Satellite State
dEvec13 Satellite::ode(const double t, const dEvec13& ode_e_state) {

	// --- PRECOMPUTATION ---

	// reinterpret const Eigen 13-vector as const State struct
	const State& ode_state = e_to_s(ode_e_state);

	// prepare derivative vector to be filled
	State derivative_of;

	// get B in Body frame
	vec3 mag_body = getMagFieldInBodyFrame(t, ode_state);

	//get v in Body frame
	vec3 v_body = glm::rotate(glm::conjugate(ode_state.q), ode_state.v);

	// compute angular-velocity-opposing magnetorque in body
	// see KPS research paper for more info
	vec3 mag_torque_body = glm::cross(-mag_gain * glm::cross(mag_body, ode_state.w), mag_body);

	double r_mag = glm::length(ode_state.r);

	vec3 aer_force_body, aer_torque_body;
	// launch aerodynamics engine
	aero.aer(aer_force_body, aer_torque_body, get1976Density(r_mag - earth.R), v_body);
	vec3 aer_force_eci = glm::rotate(ode_state.q, aer_force_body);

	// store aeroforce magnitude for later B* computation by Output engine
	aer_force_mag = glm::length(aer_force_body);

	// compute gravity gradient torque
	// see KPS reseach paper for more info
	vec3 r_body = glm::rotate(glm::conjugate(ode_state.q), ode_state.r);
	vec3 r_body_norm = glm::normalize(r_body);
	double ggt_factor = 3.0 * Earth::GM / (r_mag * r_mag * r_mag);
	// --- /PRECOMPUTATION ---


	// derivative of position is just velocity
	derivative_of.r = ode_state.v;

	// derivative of velocity is acceleration,
	// which is gravitational acceleration + aeroforce/mass
	earth.grav.getGrav(derivative_of.v, t, ode_state.r);
	derivative_of.v += aer_force_eci / m;

	// derivative of angular velocity is I^-1(tau - omega x Iomega)
	// see KPS research paper Equation 44
	// note that glm considers vec3 to be a ROW vector, so to perform
	// a matrix multiplication (I^-1) * vec, where vec should be a column
	// vector, one must actually write vec * (I^-1)
	//
	// the total torque, tau, is the sum of aerotorque, magnetorque, and
	// gravity gradient torque
	derivative_of.w = (aer_torque_body + mag_torque_body + vec3{
		ggt_factor * r_body_norm.y * r_body_norm.z * (MOI[2][2] - MOI[1][1]),
		ggt_factor * r_body_norm.z * r_body_norm.x * (MOI[0][0] - MOI[2][2]),
		ggt_factor * r_body_norm.x * r_body_norm.y * (MOI[1][1] - MOI[0][0]),
	} - glm::cross(ode_state.w, ode_state.w*MOI))*inv_MOI;

	// derivative of attitude quaternion is 1/2 (q * omega_b)
	// see KPS research paper Equation 35
	derivative_of.q = 0.5 * (ode_state.q * quat(0.0, ode_state.w));

	// reinterpret derivative State vector as an Eigen 13-vector
	return s_to_e(derivative_of);
}

// at Satellite State 'ode_state' and time 't' in seconds, obtain the B-field in the Body frame
vec3 Satellite::getMagFieldInBodyFrame(const double t, const State& ode_state) const {
	// Position: ECI -> ECEF (note: using simplified assumption of constant Earth
	// angular velocity. I use the SOFA routines for projects where extreme
	// accuracy in Earth-pointing is required.
	// I would here too, but that's simply not necessary for KPS.)
	vec3 r_ECEF = earth.eci2ecef(time_since_epoch_at_deploy + t, ode_state.r);

	// Position: ECEF -> WGS84
	double lat, lon, h;
	earth.geo.Reverse(r_ECEF.x, r_ECEF.y, r_ECEF.z, lat, lon, h);

	// Mag Field: get field at that position in ENU frame, in nanoTesla
	vec3 mag_ENU;
	earth.mag(earth.mag_yr, lat, lon, h, mag_ENU.x, mag_ENU.y, mag_ENU.z);

	// Mag Field: convert to Tesla
	mag_ENU *= TESLA_PER_nT;

	// Mag Field: (rotation only) ENU -> ECEF
	vec3 mag_ECEF = earth.enu2ecefv(lat, lon, mag_ENU);

	// Mag Field: ECEF -> ECI
	vec3 mag_ECI = earth.ecef2eci(time_since_epoch_at_deploy + t, mag_ECEF);

	//  Mag Field: ECI -> Body
	return glm::rotate(glm::conjugate(ode_state.q), mag_ECI);
}
