/*******************************************************************
*   Output.cpp
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
// Realtime or cached, ASCII or binary output of satellite
// parameters to several files, one for each parameter.
//

#include "Output.h"

// names of the output files (without file extensions, which are .bin for binary and .csv for ASCII)
const std::string Output::OUT_PREFIXES[NUM_OUTFILES]{ "t", "r", "v", "v_body", "p_e", "q", "w", "b_star", "alt"};

// name of a file created by KPS on initialization and deleted on shutdown
// (to prevent multiple instances from running simultaneously)
const std::string Output::LOCK_FILE_NAME = "KPS.lock";

bool Output::checkForExistingInstance() {
	std::ifstream lock_test(LOCK_FILE_NAME);

	if (lock_test.is_open()) {
		lock_test.close();
		std::cout << "ERROR: another instance of KPS is running! If you're sure this isn't the case," << std::endl
			<< "try again, restart your PC, or try deleting the file \"" << LOCK_FILE_NAME << '\"' << std::endl
			<< "from the KPS executable's directory. Aborting." << std::endl;
		return true;
	}

	lock_test.close();
	return false;
}

bool ASCIIOutput::init() {
	if (checkForExistingInstance()) return false;
	lock_file.open(LOCK_FILE_NAME);
	for (size_t i = 0; i < NUM_OUTFILES; ++i) {
		out[i].open(OUT_NAMES[i]);
		if (!out[i].is_open()) {
			std::cerr << "ERROR: Failed to open outfile " << OUT_NAMES[i] << ".csv! Aborting." << std::endl;
			return false;
		}
		out[i].unsetf(std::ios::floatfield);
		out[i].precision(MAX_DIGITS);
	}

	// increase performance by not constantly updating stdio ptrs
	std::ios_base::sync_with_stdio(false);
	initialized = true;
	return true;
}

bool BinaryOutput::init() {
	if (checkForExistingInstance()) return false;
	lock_file.open(LOCK_FILE_NAME);
	for (size_t i = 0; i < NUM_OUTFILES; ++i) {
		out[i].open(OUT_NAMES[i], std::ios::binary);
		if (!out[i].is_open()) {
			std::cerr << "ERROR: Failed to open outfile " << OUT_NAMES[i] << ".bin! Aborting." << std::endl;
			return false;
		}
	}

	// increase performance by not constantly updating stdio ptrs
	std::ios_base::sync_with_stdio(false);
	initialized = true;
	return true;
}

void ASCIIOutput::write(const double t, const Satellite& sat) {
	out[0] << t << '\n';
	out[1] << sat.state.r.x << ',' << sat.state.r.y << ',' << sat.state.r.z << '\n';
	out[2] << sat.state.v.x << ',' << sat.state.v.y << ',' << sat.state.v.z << '\n';

	// rotate velocity vector into Body frame
	vec3 v_body = glm::rotate(glm::conjugate(sat.state.q), sat.state.v);
	out[3] << v_body.x << ',' << v_body.y << ',' << v_body.z << '\n';

	// pointing error in degrees
	// calculated as the angle between the +z body axis and the velocity in the body frame,
	// such that if the satellite's +z is pointing in the direction of travel, error is 0
	double p_e = rad2deg(acos(v_body.z / glm::length(v_body)));
	out[4] << p_e << '\n';

	out[5] << sat.state.q.w << ',' << sat.state.q.x << ',' << sat.state.q.y << ',' << sat.state.q.z << '\n';

	out[6] << sat.state.w.x << ',' << sat.state.w.y << ',' << sat.state.w.z << '\n';

	// starred ballistic coefficient
	// see Equation 64 in KPS research paper
	double b_star = sat.aer_force_mag / (sat.m * glm::length2(sat.state.v));
	out[7] << b_star << '\n';

	// altitude in meters is distance from center of Earth minus Earth's radius (approximately)
	// a WGS84 geoid height could be used instead for more accuracy.
	// assuming spherical isolates orbital behavior in the plot, which is desirable for troubleshooting
	double alt = glm::length(sat.state.r) - Earth::R;
	out[7] << alt << '\n';

	// the OS will cache writes unless flushed,
	// which decreases performance but is desirable for realtime plotting
	if (realtime) {
		for (int i = 0; i < NUM_OUTFILES; ++i) out[i].flush();
	}
}

void BinaryOutput::write(const double t, const Satellite& sat) {
	out[0].write(reinterpret_cast<const char*>(&t), sizeof(double));
	out[1].write(reinterpret_cast<const char*>(&sat.state.r), sizeof(vec3));
	out[2].write(reinterpret_cast<const char*>(&sat.state.v), sizeof(vec3));

	// rotate velocity vector into Body frame
	vec3 v_body = glm::rotate(glm::conjugate(sat.state.q), sat.state.v);
	out[3].write(reinterpret_cast<const char*>(&v_body), sizeof(vec3));

	// pointing error in degrees
	// calculated as the angle between the +z body axis and the velocity in the body frame,
	// such that if the satellite's +z is pointing in the direction of travel, error is 0
	double p_e = rad2deg(acos(v_body.z / glm::length(v_body)));
	out[4].write(reinterpret_cast<const char*>(&p_e), sizeof(double));

	out[5].write(reinterpret_cast<const char*>(&sat.state.q), sizeof(quat));

	out[6].write(reinterpret_cast<const char*>(&sat.state.w), sizeof(vec3));

	// starred ballistic coefficient
	// see Equation 64 in KPS research paper
	double b_star = sat.aer_force_mag / (sat.m * glm::length2(sat.state.v));
	out[7].write(reinterpret_cast<const char*>(&b_star), sizeof(double));

	// altitude in meters is distance from center of Earth minus Earth's radius (approximately)
	// a WGS84 geoid height could be used instead for more accuracy.
	// assuming spherical isolates orbital behavior in the plot, which is desirable for troubleshooting
	double alt = glm::length(sat.state.r) - Earth::R;
	out[8].write(reinterpret_cast<const char*>(&alt), sizeof(double));

	// the OS will cache writes unless flushed,
	// which decreases performance but is desirable for realtime plotting
	if (realtime) {
		for (int i = 0; i < NUM_OUTFILES; ++i) out[i].flush();
	}
}

Output::~Output() {
	lock_file.close();
	std::remove(LOCK_FILE_NAME.c_str());
	if (initialized) {
		for (size_t i = 0; i < NUM_OUTFILES; ++i) out[i].close();
	}
}