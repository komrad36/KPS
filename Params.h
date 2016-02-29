/*******************************************************************
*   Params.h
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
// Interprets map of raw parameter name and value strings into
// valid parameters, including performing validity checks.
//

#pragma once

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

#include "Earth.h"
#include "Aero.h"
#include "Propagators.h"
#include "Output.h"
#include "glm_util.h"
#include "Eigen_util.h"
#include "str_util.h"

typedef std::unordered_map<std::string, std::string> key_val_map;
typedef std::pair<std::string, std::string> key_val_pair;

class Params {
public:

	static const std::vector<std::string> param_names;
	static const int NUM_PARAMS = 21;

	static const int MOI_ELEMENTS = 9;
	static const int SAT_CM_ELEMENTS = 3;
	static const int SAT_INIT_POS_ELEMENTS = 3;
	static const int SAT_INIT_Q_ELEMENTS = 4;
	static const int SAT_INIT_V_ELEMENTS = 3;
	static const int SAT_INIT_W_ELEMENTS = 3;
	
	double mag_gain;
	int aero_mode;
	double time_since_epoch_at_deploy;
	std::string grav_model;
	Propagator<ODE_VEC_N>* propagator;
	std::string propagator_name;
	double abs_tol;
	double rel_tol;
	double max_step_size;
	std::string mag_model;
	double mag_year;
	double pitch;
	Output* output;
	vec3 sat_cm;
	vec3 sat_init_pos;
	quat sat_init_q;
	vec3 sat_init_v;
	vec3 sat_init_w;
	double sat_mass;
	mat3 sat_moi;
	double time_span;
	bool binary_output;
	bool realtime_output;

private:

public:
	Params() : propagator(nullptr), output(nullptr) {}

	~Params() {
		if (propagator != nullptr) delete propagator;
		if (output != nullptr) delete output;
	}

	bool assign(const key_val_map& raw_params);
	void printVals();

	static void printNames();

private:

	static const std::vector<std::string> grav_names;
	static const std::vector<std::string> mag_names;

};


