/*******************************************************************
*   main.cpp
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
// KPS is a simultaneous orbital and attitude propagator for
// satellites in Low-Earth Orbit, using either CUDA or CPU
// for aerodynamics simulation. See my GitHub, the README file,
// and/or the KPS research paper for more information.
//


// milliseconds between prints of current propagator step and time
// < 40 might result in performance impact, especially on Windows
const double MS_PER_PRINT = 40.0;

#define KPS_VER "1.0"

const double SEC_PER_MS = 0.001;
const int EXPECTED_ARGS = 2;

// signal handling to gracefully deal with CTRL+C or console exit
#ifdef _WIN32
#include <Windows.h>
#define ADD_HANDLER static_cast<BOOL>(1)
#else
#include <signal.h>
#include <unistd.h>
#endif

#include <chrono>
#include <functional>
#include <iostream>
#include <vector>

#include "Earth.h"
#include "KAero.h"
#include "ParamInput.h"
#include "PolyInput.h"
#include "Output.h"
#include "Params.h"
#include "Propagators.h"
#include "Satellite.h"

const std::string CUDA_EULA_FILE_NAME = "CUDA_EULA_ACCEPTED";

using std::chrono::system_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

volatile bool is_closing;

// signal handling to gracefully deal with CTRL+C or console exit
#ifdef _WIN32
BOOL CtrlHandler(DWORD fdwCtrlType) {
	// silence unused parameter
	static_cast<void>(fdwCtrlType);

	// tell main thread to close
	is_closing = true;

	// don't let Windows kill the process
	Sleep(INFINITE);

	// will never be reached; silence compiler about no return path
	return true;
}
#else
void CtrlHandler(int s) {
	// silence unused parameter
	static_cast<void>(s);	

	// Linux handles things differently.
	// as long as you've registered the handler,
	// it won't kill the process,
	// so just tell main thread to close
	is_closing = true;
}
#endif

// if debugging in VS on Windows, check for memory leaks
#if defined(_DEBUG) && defined(_WIN32)

#ifndef DBG_NEW
#define DBG_NEW new ( _NORMAL_BLOCK , __FILE__ , __LINE__ )
#define new DBG_NEW
#endif // DBG_NEW

#ifndef _CRTDBG_MAP_ALLOC
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif // _CRTDBG_MAP_ALLOC

#endif // _DEBUG

void printHeader() {
	std::cout
		<< "# ----------------------- #\n"
		<< "#           KPS           #\n"
		<< "# ----------------------- #\n"
		<< "#       Version " << KPS_VER << "       #\n"
		<< "# ----------------------- #\n"
		<< "#       Kareem Omar       #\n"
		<< "#   kareem.omar@uah.edu   #\n"
		<< "#   github.com/komrad36   #\n"
		<< "# ----------------------- #\n"
		<< std::endl;
}

void printUsage() {
	std::cout
		<< "Usage: KPS <polygon_file> <parameters_file>\n"
		<< "where <polygon_file> is the name of\n"
		<< "a polygon input file containing comma-separated\n"
		<< "vertices, one per line, with every " << NUM_VTX << "\n"
		<< "vertices making up one polygon describing a panel\n"
		<< "of the satellite.\n\n"
		<< "<parameters_file> is a KPS input file containing\n"
		<< "name-value pairs, in any order, one per line, such as:\n"
		<< "MAG_GAIN = 15000\n"
		<< "Spaces and capitalization are ignored.\n"
		<< "See the README for more.\n"
		<< "Required parameters are:\n"
		<< std::endl;
	Params::printNames();
}

int main(int argc, char* argv[]) {

	is_closing = false;

	// set FPU to FTZ and DAZ to eliminate
	// performance issues due to denormal handling
	// see https://en.wikipedia.org/wiki/Denormal_number
	_MM_SET_FLUSH_ZERO_MODE(_MM_FLUSH_ZERO_ON);
	_MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);

	// if debugging in VS on Windows, check for memory leaks
#if defined(_DEBUG) && defined(_WIN32)
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif

	printHeader();

	if (argc != EXPECTED_ARGS + 1) {
		printUsage();
		return EXIT_FAILURE;
	}

	// --- Parse polygon file ---
	std::cout
		<< "KPS Initializing..." << std::endl
		<< "Parsing requested polygon file \"" << argv[1] << "\"... ";
	// satellite geometry is an array of 3-vectors,
	// with each 3-vector representing a vertex
	std::vector<vec3> poly;
	int num_poly;
	if (!parsePoly(poly, num_poly, argv[1])) return EXIT_FAILURE;
	std::cout << num_poly << " polygons loaded." << std::endl;
	// --- /Parse polygon file ---


	// --- Parse parameters file ---
	// a map of param names to param values
	key_val_map raw_params;
	std::cout << "Parsing requested parameters file \"" << argv[2] << "\"... ";
	if (!parseParams(raw_params, argv[2])) return EXIT_FAILURE;
	Params params;
	if (!params.assign(raw_params)) return EXIT_FAILURE;
	std::cout
		<< params.NUM_PARAMS << " parameters loaded." << std::endl << std::endl
		<< "-------------------------------- Configuration --------------------------------" << std::endl;
	params.printVals();
	std::cout << std::string(79, '-') << std::endl << std::endl;
	// --- Parse parameters file ---


	// --- Initialize gravity model ---
	std::cout << "Initializing gravitational model..." << std::endl;
	GravityModel* grav_model;
	// GeographicLib constructors throw a GeographicException on failure to init
	try {
		if (params.grav_model == "point") {
			grav_model = new PointGravityModel();
		}
		else if (params.grav_model == "wgs84") {
			grav_model = new GeographicLibWGS84GravityModel();
		}
		else {
			grav_model = new GeographicLibOtherGravityModel(params.grav_model);
		}
	}
	catch (const std::exception& e) {
		std::cerr << "ERROR: " << e.what() << std::endl
			<< "Please request a different gravity model or download the requested one from" << std::endl
			<< "http://geographiclib.sourceforge.net/html/gravity.html#gravityinst" << std::endl;
		return EXIT_FAILURE;
	}
	// --- /Initialize gravity model ---


	// --- Initialize magnetic model ---
	std::cout << "Initializing magnetic model..." << std::endl;
	GeographicLib::MagneticModel* mag_model;
	// GeographicLib constructors throw a GeographicException on failure to init
	try {
		mag_model = new GeographicLib::MagneticModel{ params.mag_model };
	}
	catch (const std::exception& e) {
		std::cerr << "ERROR: " << e.what() << std::endl
			<< "Please request a different magnetic model or download the requested one from" << std::endl
			<< "http://geographiclib.sourceforge.net/html/magnetic.html#magneticinst" << std::endl;
		return EXIT_FAILURE;
	}
	// --- /Initialize magnetic model ---


	// --- Initialize geoid model ---
	std::cout << "Initializing geoid model..." << std::endl;
	const Earth earth(params.mag_year, *grav_model, *mag_model);
	// --- /Initialize geoid model ---

	// --- Initialize aero engine ---
	std::cout << "Initializing KAero..." << std::endl;
	KAero* kdrag;
	if (params.cuda_device == USE_CPU) {
		kdrag = new KAero_CPU(params.pitch, num_poly, &poly[0], params.sat_cm);
		std::cout << "CPU KAero ready." << std::endl;
	}
	else {
		// if in CUDA mode, check if the user has accepted the CUDA EULA
		std::ifstream eula_check(CUDA_EULA_FILE_NAME);
		if (!eula_check.is_open()) {
			// if not, prompt the user immediately.
			std::cout
				<< "\nKPS has been started in CUDA mode. CUDA is (c) NVIDIA Corporation.\n"
				<< "Use of the CUDA runtime libraries requires that you agree to the CUDA End-User\n"
				<< "License Agreement (EULA) available at http://docs.nvidia.com/cuda/eula/\n\n"
				<< "Do you agree to these terms? Type YES to indicate agreement.\n" << std::endl;

			std::string response;
			std::cin >> response;
			std::cout << std::endl;
			
			if (response == "YES") {
				// if user accepts EULA, record user's choice for future runs
				std::ofstream eula_save(CUDA_EULA_FILE_NAME);
				eula_save.close();
			}
			else {
				// if not, terminate.
				std::cerr
					<< "ERROR: user did not accept the NVIDIA CUDA EULA. KPS can only be run in CPU\n"
					<< "mode unless the user accepts the NVIDIA CUDA EULA. Aborting." << std::endl;
				return EXIT_FAILURE;
			}

		}
		eula_check.close();

		KAero_CUDA* kdrag_cuda = new KAero_CUDA(params.pitch, num_poly, &poly[0], params.sat_cm);
		if (!kdrag_cuda->init(params.cuda_device)) return EXIT_FAILURE;
		kdrag = kdrag_cuda;
		std::cout << "CUDA KAero ready." << std::endl;
	}
	// polygon data has been copied into aero engine and is no longer needed
	poly.clear();
	// --- /Initialize aero engine ---


	// --- Set signal handler on either OS ---
#ifdef _WIN32
	SetConsoleCtrlHandler(reinterpret_cast<PHANDLER_ROUTINE>(CtrlHandler), ADD_HANDLER);
#else
	struct sigaction sigIntHandler;
	sigIntHandler.sa_handler = CtrlHandler;
	sigemptyset(&sigIntHandler.sa_mask);
	sigIntHandler.sa_flags = 0;
	sigaction(SIGINT, &sigIntHandler, NULL);
	sigaction(SIGTERM, &sigIntHandler, NULL);
	sigaction(SIGHUP, &sigIntHandler, NULL);
#endif
	// --- /Set signal handler on either OS ---


	// --- Initialize output engine ---
	std::cout << "Initializing output engine..." << std::endl;
	Output& output = *params.output;
	if (!output.init()) return EXIT_FAILURE;
	// --- /Initialize output engine ---


	// --- Initialize satellite model ---
	std::cout << "Spawning satellite..." << std::endl;
	Satellite sat(earth, *kdrag, params.sat_mass, params.mag_gain, params.time_since_epoch_at_deploy,
	{ params.sat_init_pos, params.sat_init_v, params.sat_init_q, params.sat_init_w },
	params.sat_moi);
	// --- /Initialize satellite model ---


	// --- Initialize propagator ---
	std::cout << "Initializing propagator..." << std::endl;
	Propagator<ODE_VEC_N>& propagator(*params.propagator);
	double t = 0.0;
	propagator.init(std::bind(&Satellite::ode, sat, std::placeholders::_1, std::placeholders::_2), t, sat.e_state, params.abs_tol, params.rel_tol, params.max_step_size);
	output.write(t, sat);
	// --- /Initialize propagator ---

	std::cout << std::endl
		<< "KPS READY." << std::endl << std::endl
		<< "Propagating:" << std::endl;

	// record the start time and prepare to track time between prints
	system_clock::time_point current_RTC, last_print_RTC, start_RTC = system_clock::now();

	// no platform-independent way to printf() a size_t, so using unsigned long long instead
	unsigned long long steps = 0;

	// --- Main Propagation Loop ---
	while (propagator.step(t, sat.e_state) < params.time_span && sat.isInValidState()) {

		// optional renormalization of quaternion to minimize the impact of integrator numerical error
		// not that crucial because the integrators are really good
		// in one test, the quaternion norm was 1.00000000596 after 2 years of simulated orbit (!)
		sat.state.q = glm::normalize(sat.state.q);

		output.write(t, sat);
		++steps;
		if (duration_cast<milliseconds>((current_RTC = system_clock::now()) - last_print_RTC).count() > MS_PER_PRINT) {
			last_print_RTC = current_RTC;
			printf("\r Step %llu | %f sec       ", steps, t);
			fflush(stdout);
			if (is_closing) break;
		}
	}
	// --- /Main Propagation Loop ---

	if (!is_closing) printf("\r Step %llu | %f sec         \n", ++steps, sat.isInValidState() ? params.time_span : t);
	fflush(stdout);

	std::cout
		<< (is_closing ? "\nPropagation terminated by user after " :
		sat.isInValidState() ? "Propagation complete in " : "Satellite DEORBITED! Propagation complete in ")
		<< static_cast<double>(duration_cast<milliseconds>(system_clock::now() - start_RTC).count()) * SEC_PER_MS
		<< " sec realtime." << std::endl << std::endl;

	delete kdrag;
	delete grav_model;
	delete mag_model;

}