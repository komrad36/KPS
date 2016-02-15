/*******************************************************************
*   Output.h
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

#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream>
#include <iostream>

#include "Earth.h"
#include "Satellite.h"

// pure virtual base class to allow branchless selection
// of binary or ASCII output
class Output {
protected:
	bool realtime;

	static const int NUM_OUTFILES = 9;
	static const std::string OUT_PREFIXES[NUM_OUTFILES];
	static const std::string LOCK_FILE_NAME;

	std::ofstream out[NUM_OUTFILES];

	bool initialized;

	bool checkForExistingInstance();

	std::ofstream lock_file;

public:
	virtual void write(const double t, const Satellite& sat) = 0;
	virtual bool init() = 0;

	Output(bool realtime_output) : realtime(realtime_output), initialized(false) {}
	virtual ~Output() = 0;
};

class ASCIIOutput : public Output {
private:
	std::string OUT_NAMES[NUM_OUTFILES];
public:
	ASCIIOutput(bool realtime_output) : Output(realtime_output) {
		for (int i = 0; i < NUM_OUTFILES; ++i) {
			OUT_NAMES[i] = OUT_PREFIXES[i] + ".csv";
		}
	}

	bool init();
	void write(const double t, const Satellite& sat);
};

class BinaryOutput : public Output {
private:
	std::string OUT_NAMES[NUM_OUTFILES];
public:
	BinaryOutput(bool realtime_output) : Output(realtime_output) {
		for (int i = 0; i < NUM_OUTFILES; ++i) {
			OUT_NAMES[i] = OUT_PREFIXES[i] + ".bin";
		}
	}

	bool init();
	void write(const double t, const Satellite& sat);
};

#endif
