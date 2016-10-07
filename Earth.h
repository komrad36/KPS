/*******************************************************************
*   Earth.h
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
// Earth (orbited body) model, including wrappers for gravitational modeling
//

#pragma once

#include <GeographicLib/Geocentric.hpp>
#include <GeographicLib/GravityModel.hpp>
#include <GeographicLib/MagneticModel.hpp>

#include <string>

#include "glm_util.h"
#include "Eigen_util.h"

// pure virtual base class for gravity wrappers
struct GravityModel {
	virtual void getGrav(vec3& ret, const double t, const vec3& pos) const = 0;
	virtual ~GravityModel() {}
};

class Earth {
public:

	// rate of Earth's rotation about its axis [rad/sec]
	static constexpr double Earth::OMEGA = 7.29211509e-5;

	// standard gravitational parameter
	static constexpr double Earth::GM = 398600441800000.0;

	// Earth's radius [m]
	static constexpr double Earth::R = 6371000.0;

	// altitude at which deorbit should be assumed [m]
	static constexpr double Earth::DEORBIT_ALT = 100000.0;

	// Earth's magnetic field is volatile, so the desired year is required as an input to GeographicLib
	// magnetic models
	const double mag_yr;

	const GravityModel& grav;
	const GeographicLib::MagneticModel& mag;
	const GeographicLib::Geocentric& geo;

	// convert a vector in the ECI frame to the ECEF frame
	// the two frames are identical when t = 0
	static inline vec3 eci2ecef(const double t_since_epoch, const vec3& eci) {
		double yaw = t_since_epoch * OMEGA;
		double sin_yaw = sin(yaw), cos_yaw = cos(yaw);
		return { cos_yaw*eci.x + sin_yaw*eci.y, -sin_yaw*eci.x + cos_yaw*eci.y, eci.z };
	}

	// convert a vector in the ECEF frame to the ECI frame
	// the two frames are identical when t = 0
	static inline vec3 ecef2eci(const double t_since_epoch, const vec3& ecef) {
		double yaw = t_since_epoch * OMEGA;
		double sin_yaw = sin(yaw), cos_yaw = cos(yaw);
		return { cos_yaw*ecef.x - sin_yaw*ecef.y, sin_yaw*ecef.x + cos_yaw*ecef.y, ecef.z };
	}

	// ROTATE a vector from the ENU frame into the ECEF frame
	// NOTE: this is NOT the same as converting an ENU position into an ECEF position;
	// only the orientation is considered in this function, as required to convert
	// velocities or magnetic field components, for example.
	static inline vec3 enu2ecefv(const double lat_deg, const double lon_deg, const vec3& enu) {
		double lat = deg2rad(lat_deg);
		double lon = deg2rad(lon_deg);
		double sin_lat = sin(lat);
		double cos_lat = cos(lat);
		double sin_long = sin(lon);
		double cos_long = cos(lon);
		double precomp = cos_lat*enu.z - sin_lat*enu.y;
		return{
			cos_long*precomp - sin_long*enu.x,
			sin_long*precomp + cos_long*enu.x,
			cos_lat*enu.y + sin_lat*enu.z
		};
	}

public:
	Earth(double mag_year, const GravityModel& gravity_model, const GeographicLib::MagneticModel& magnetic_model) :
		mag_yr(mag_year), grav(gravity_model), mag(magnetic_model), geo{ GeographicLib::Geocentric::WGS84() } {}

};

// point mass (or equivalently, uniform spherical Earth) gravity derived class
struct PointGravityModel : public GravityModel {
	PointGravityModel() {}

	inline void getGrav(vec3& ret, const double t, const vec3& pos) const {
		// silence unused parameter
		static_cast<void>(t);

		// a = GM/|r|^2 in the direction of -r
		ret = -Earth::GM / glm::length2(pos) * glm::normalize(pos);
	}
};

// WGS84 gravity
// This particular model gets its own class because it's the only
// model that is rotationally symmetric about the axis of rotation
// (the ECI and ECEF z-axis), which means it is NOT necessary
// to first convert the ECI position to ECEF, then get
// gravity, then convert back, as it is with other models,
// allowing an increase in performance without branching.
struct GeographicLibWGS84GravityModel : public GravityModel {
	const GeographicLib::GravityModel grav;

	GeographicLibWGS84GravityModel() : grav{ "wgs84" } {}
	inline void getGrav(vec3& ret, const double t, const vec3& pos) const {
		static_cast<void>(t);
		grav.V(pos.x, pos.y, pos.z, ret.x, ret.y, ret.z);
	}
};

// any other gravity model such as EGM96
struct GeographicLibOtherGravityModel : public GravityModel {
	const GeographicLib::GravityModel grav;

	GeographicLibOtherGravityModel(const std::string& gravity_model) : grav{gravity_model} {}

	inline void getGrav(vec3& ret, const double t, const vec3& pos) const {

		// position ECI -> ECEF
		vec3 pos_ecef = Earth::eci2ecef(t, pos);

		// get gravity in ECEF
		double grav_eci_x, grav_eci_y, grav_eci_z;
		grav.V(pos_ecef.x, pos_ecef.y, pos_ecef.z, grav_eci_x, grav_eci_y, grav_eci_z);

		// gravity ECEF -> ECI
		ret = Earth::ecef2eci(t, {grav_eci_x, grav_eci_y, grav_eci_z});
	}
};


