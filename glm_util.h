/*******************************************************************
*   glm_util.h
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
// Includes, constants, and typedefs for KPS modules that deal with
// vectors, polygons, and matrices using the GLM library
//

#pragma once

#define NUM_VTX (4)

#include <limits>

#include <glm/glm.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/norm.hpp>

// number of elements in the vectors
// integrated by the ODE solvers
//
// for orbital and attitude propagation, it's 13:
// 3 for position
// 3 for velocity
// 4 for orientation (quaternion)
// 3 for angular velocity
constexpr size_t ODE_VEC_N = 13;

constexpr int AUTO_SELECT = -1;

constexpr int GRID = -2;

constexpr int ANALYTICAL = -3;

constexpr double TESLA_PER_nT = 1e-9;

constexpr size_t STD_DIGITS = std::numeric_limits<double>::digits10;

constexpr size_t MAX_DIGITS = std::numeric_limits<double>::max_digits10;

// KPS deals in very small numbers with long accumulations,
// especially in the integrators.
// Single-precision
// is NOT sufficient.
typedef glm::dmat3 mat3;
typedef glm::dvec3 vec3;
typedef glm::dquat quat;

const vec3 X_HAT{ 1.0, 0.0, 0.0 };
const vec3 Y_HAT{ 0.0, 1.0, 0.0 };
const vec3 Z_HAT{ 0.0, 0.0, 1.0 };

// deg * PI / 180
constexpr inline double deg2rad(const double deg) { return deg * 0.0174532925199432957; }

// rad * 180 / PI
constexpr inline double rad2deg(const double rad) { return rad * 57.295779513082320876; }

