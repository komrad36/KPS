/*******************************************************************
*   Params.cpp
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
// Interprets map of raw parameter name and value strings into
// valid parameters, including performing validity checks.
//

#include "Params.h"

// names of supported gravitational models
const std::vector<std::string> Params::grav_names{ "POINT", "WGS84", "EGM84", "EGM96", "EGM2008" };

// names of supported magnetic models
const std::vector<std::string> Params::mag_names{ "WMM2010", "WMM2015", "IGRF11", "IGRF12", "EMM2010", "EMM2015" };

// names of all required parameters
const std::vector<std::string> Params::param_names{ "MAG_GAIN", "CUDA_DEVICE", "TIME_SINCE_EPOCH_AT_DEPLOY", "GRAV_MODEL", "PROPAGATOR",
													"ABS_TOL", "REL_TOL", "KAERO_PITCH", "MAX_STEP_SIZE", "MAG_MODEL", "MAG_YEAR",
													"BINARY_OUTPUT", "REALTIME_OUTPUT", "SAT_CM", "SAT_INIT_POS", "SAT_INIT_Q",
													"SAT_INIT_V", "SAT_INIT_W", "SAT_MASS", "SAT_MOI", "TIME_SPAN" };

// assign specific values to varaibles in parameter struct based on raw strings from ParamInput
bool Params::assign(const key_val_map& raw_params) {
	key_val_map::const_iterator pair;
	const std::string err = "ERROR: One or more required parameters missing. Aborting.\nSee documentation for required parameters.";

	// --- Magnetorque Gain ---
	if ((pair = raw_params.find("MAG_GAIN")) != raw_params.end()) {
		mag_gain = atof(pair->second.c_str());
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Magnetorque Gain ---


	// --- Satellite Center of Mass ---
	if ((pair = raw_params.find("SAT_CM")) != raw_params.end()) {
		std::vector<std::string> vals;
		tokenize(vals, pair->second, ',');
		if (vals.size() != SAT_CM_ELEMENTS) {
			std::cerr << std::endl << "ERROR: center of mass requires " << SAT_CM_ELEMENTS << " elements. Found " << vals.size() << ". Aborting." << std::endl;
			return false;
		}
		else {
			for (int i = 0; i < SAT_CM_ELEMENTS; ++i) sat_cm[i] = atof(vals[i].c_str());
		}
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Satellite Center of Mass ---


	// --- Aerodynamics Simulation Linear Pitch ---
	if ((pair = raw_params.find("KAERO_PITCH")) != raw_params.end()) {
		pitch = atof(pair->second.c_str());
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Aerodynamics Simulation Linear Pitch ---


	// --- CUDA Device ID ---
	if ((pair = raw_params.find("CUDA_DEVICE")) != raw_params.end()) {
		if (pair->second == "AUTO") {
			cuda_device = AUTO_SELECT;
		}
		else if (pair->second == "NONE") {
			cuda_device = USE_CPU;
		}
		else if ((cuda_device = atoi(pair->second.c_str())) < 0) {
			std::cerr << std::endl << "ERROR: for CUDA_DEVICE, specify a valid CUDA device ID," << std::endl
				<< "or AUTO to auto-select the device with the highest GFLOPS," << std::endl
				<< "or NONE to disable CUDA and use the CPU exclusively." << std::endl;
			return false;
		}
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /CUDA Device ID ---


	// --- Seconds Since Epoch at Satellite Deploy ---
	if ((pair = raw_params.find("TIME_SINCE_EPOCH_AT_DEPLOY")) != raw_params.end()) {
		time_since_epoch_at_deploy = atof(pair->second.c_str());
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Seconds Since Epoch at Satellite Deploy ---


	// --- Gravitational Model ---
	if ((pair = raw_params.find("GRAV_MODEL")) != raw_params.end()) {
		if (std::any_of(grav_names.begin(), grav_names.end(), [pair](std::string s){return pair->second == s; })) {
			grav_model = pair->second;
			tolower(grav_model);
		}
		else {
			std::cerr << std::endl << "ERROR: specified gravity model not found. Aborting." << std::endl;
			return false;
		}
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Gravitational Model ---


	// --- Propagator ---
	if ((pair = raw_params.find("PROPAGATOR")) != raw_params.end()) {
		if (pair->second == "RKDP") {
			propagator = new RKDP<ODE_VEC_N>();
		}
		else if (pair->second == "ABM") {
			propagator = new ABM<ODE_VEC_N>();
		}
		else {
			std::cerr << std::endl << "ERROR: specified propagator not found. Aborting." << std::endl;
			return false;
		}
		propagator_name = pair->second;
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Propagator ---


	// --- Integrator Absolute Tolerance ---
	if ((pair = raw_params.find("ABS_TOL")) != raw_params.end()) {
		abs_tol = atof(pair->second.c_str());
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Integrator Absolute Tolerance ---


	// --- Integrator Relative Tolerance ---
	if ((pair = raw_params.find("REL_TOL")) != raw_params.end()) {
		rel_tol = atof(pair->second.c_str());
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Integrator Absolute Tolerance ---


	// --- Integrator Maximum Step Size ---
	if ((pair = raw_params.find("MAX_STEP_SIZE")) != raw_params.end()) {
		max_step_size = atof(pair->second.c_str());
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Integrator Maximum Step Size ---


	// --- Magnetic Model ---
	if ((pair = raw_params.find("MAG_MODEL")) != raw_params.end()) {
		if (std::any_of(mag_names.begin(), mag_names.end(), [pair](std::string s){return pair->second == s; })) {
			mag_model = pair->second;
			tolower(mag_model);
		}
		else {
			std::cerr << std::endl << "ERROR: specified magnetic model not found. Aborting." << std::endl;
			return false;
		}
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Magnetic Model ---


	// --- Magnetic Model Year ---
	if ((pair = raw_params.find("MAG_YEAR")) != raw_params.end()) {
		mag_year = atof(pair->second.c_str());
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Magnetic Model Year ---


	// --- Binary Output? ---
	if ((pair = raw_params.find("BINARY_OUTPUT")) != raw_params.end()) {
		binary_output = (pair->second == "TRUE");
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Binary Output? ---


	// --- Realtime Output? ---
	if ((pair = raw_params.find("REALTIME_OUTPUT")) != raw_params.end()) {
		realtime_output = (pair->second == "TRUE");
		if (binary_output) {
			output = new BinaryOutput(realtime_output);
		}
		else {
			output = new ASCIIOutput(realtime_output);
		}
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Realtime Output? ---


	// --- Satellite Initial Position ---
	if ((pair = raw_params.find("SAT_INIT_POS")) != raw_params.end()) {
		std::vector<std::string> vals;
		tokenize(vals, pair->second, ',');
		if (vals.size() != SAT_INIT_POS_ELEMENTS) {
			std::cerr << std::endl << "ERROR: initial position requires " << SAT_INIT_POS_ELEMENTS << " elements. Found " << vals.size() << ". Aborting." << std::endl;
			return false;
		}
		else {
			for (int i = 0; i < SAT_INIT_POS_ELEMENTS; ++i) sat_init_pos[i] = atof(vals[i].c_str());
		}
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Satellite Initial Position ---


	// --- Satellite Initial Orientation ---
	if ((pair = raw_params.find("SAT_INIT_Q")) != raw_params.end()) {
		std::vector<std::string> vals;
		tokenize(vals, pair->second, ',');
		if (vals.size() != SAT_INIT_Q_ELEMENTS) {
			std::cerr << std::endl << "ERROR: initial quaternion requires " << SAT_INIT_Q_ELEMENTS << " elements. Found " << vals.size() << ". Aborting." << std::endl;
			return false;
		}
		else {
			for (int i = 0; i < SAT_INIT_Q_ELEMENTS; ++i) sat_init_q[i ? i - 1 : 3] = atof(vals[i].c_str());
		}
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Satellite Initial Orientation ---


	// --- Satellite Initial Velocity ---
	if ((pair = raw_params.find("SAT_INIT_V")) != raw_params.end()) {
		std::vector<std::string> vals;
		tokenize(vals, pair->second, ',');
		if (vals.size() != SAT_INIT_V_ELEMENTS) {
			std::cerr << std::endl << "ERROR: initial velocity requires " << SAT_INIT_V_ELEMENTS << " elements. Found " << vals.size() << ". Aborting." << std::endl;
			return false;
		}
		else {
			for (int i = 0; i < SAT_INIT_V_ELEMENTS; ++i) sat_init_v[i] = atof(vals[i].c_str());
		}
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Satellite Initial Velocity ---


	// --- Satellite Initial Angular Velocity ---
	if ((pair = raw_params.find("SAT_INIT_W")) != raw_params.end()) {
		std::vector<std::string> vals;
		tokenize(vals, pair->second, ',');
		if (vals.size() != SAT_INIT_W_ELEMENTS) {
			std::cerr << std::endl << "ERROR: initial velocity requires " << SAT_INIT_W_ELEMENTS << " elements. Found " << vals.size() << ". Aborting." << std::endl;
			return false;
		}
		else {
			for (int i = 0; i < SAT_INIT_W_ELEMENTS; ++i) sat_init_w[i] = atof(vals[i].c_str());
		}
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Satellite Initial Angular Velocity ---


	// --- Satellite Mass ---
	if ((pair = raw_params.find("SAT_MASS")) != raw_params.end()) {
		sat_mass = atof(pair->second.c_str());
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Satellite Mass ---


	// --- Satellite Moment of Inertia Matrix ---
	if ((pair = raw_params.find("SAT_MOI")) != raw_params.end()) {
		std::vector<std::string> vals;
		tokenize(vals, pair->second, ',');
		if (vals.size() != MOI_ELEMENTS) {
			std::cerr << std::endl << "ERROR: MOI matrix requires " << MOI_ELEMENTS << " elements. Found " << vals.size() << ". Aborting." << std::endl;
			return false;
		}
		else {
			for (int i = 0; i < MOI_ELEMENTS; ++i) sat_moi[i / 3][i % 3] = atof(vals[i].c_str());
		}
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}
	// --- /Satellite Moment of Inertia Matrix ---


	// --- Time Span of Simulation ---
	if ((pair = raw_params.find("TIME_SPAN")) != raw_params.end()) {
		time_span = atof(pair->second.c_str());
	}
	else {
		std::cerr << std::endl << err << std::endl;
		return false;
	}

	if (raw_params.size() > NUM_PARAMS) {
		std::cerr << std::endl << "WARN: Ignoring one or more unrecognized parameters." << std::endl;
	}
	// --- /Time Span of Simulation ---


	return true;
}

void Params::printNames() {
	for (auto&& param : param_names) {
		std::cout << param << std::endl;
	}
	std::cout << std::endl;
}

void Params::printVals() {
	std::cout.precision(STD_DIGITS);
	std::cout
		<< "Magnetorque Gain: " << mag_gain << std::endl
		<< (cuda_device == USE_CPU ? "CPU Mode" : "CUDA Mode - Device ID: " + (cuda_device == -1 ? "AUTO" : std::to_string(cuda_device))) << std::endl
		<< "Time Since Epoch at Deploy: " << time_since_epoch_at_deploy << " s" << std::endl
		<< "Gravitational Model: " << grav_model << std::endl
		<< "Propagator: " << propagator_name << std::endl
		<< "Absolute Tolerance: " << abs_tol << std::endl
		<< "Relative Tolerance: " << rel_tol << std::endl
		<< "Max Step Size: " << max_step_size << " s" << std::endl
		<< "KAero Linear Pitch: " << pitch << " m" << std::endl
		<< "Magnetic Model: " << mag_model << std::endl
		<< "Magnetic Model Year: " << mag_year << std::endl
		<< (binary_output ? "Binary" : "ASCII") << " Output" << std::endl
		<< (realtime_output ? "Realtime" : "Non-realtime") << " Output" << std::endl
		<< "Sat Initial Position: { " << sat_init_pos[0] << ", " << sat_init_pos[1] << ", " << sat_init_pos[2] << " }" << " m" << std::endl
		<< "Sat Initial Attitude: { " << sat_init_q.w << ", " << sat_init_q.x << ", " << sat_init_q.y << ", " << sat_init_q.z << " }" << std::endl
		<< "Sat Initial Velocity: { " << sat_init_v[0] << ", " << sat_init_v[1] << ", " << sat_init_v[2] << " }" << " m/s" << std::endl
		<< "Sat Initial Ang Vel:  { " << sat_init_w[0] << ", " << sat_init_w[1] << ", " << sat_init_w[2] << " }" << " rad/s" << std::endl
		<< "Sat Center of Mass:   { " << sat_cm[0] << ", " << sat_cm[1] << ", " << sat_cm[2] << " }" << " m" << std::endl
		<< "Sat Mass: " << sat_mass << " kg" << std::endl
		<< "Sat MOI: " << std::setw(MAX_DIGITS) << sat_moi[0][0] << ", " << std::setw(MAX_DIGITS) << sat_moi[0][1] << ", " << std::setw(MAX_DIGITS) << sat_moi[0][2] << std::endl
		<< "         " << std::setw(MAX_DIGITS) << sat_moi[1][0] << ", " << std::setw(MAX_DIGITS) << sat_moi[1][1] << ", " << std::setw(MAX_DIGITS) << sat_moi[1][2] << std::endl
		<< "         " << std::setw(MAX_DIGITS) << sat_moi[2][0] << ", " << std::setw(MAX_DIGITS) << sat_moi[2][1] << ", " << std::setw(MAX_DIGITS) << sat_moi[2][2] << std::endl
		<< "Time Span: " << time_span << " s" << std::endl;
}
