/*******************************************************************
*   ParamInput.cpp
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
// Parses a KPS parameters file into a map of parameter names and values.
//

#include "ParamInput.h"

bool parseParams(key_val_map& params, const std::string& filename, const char delim_char, const char comment_char) {
	std::ifstream in_f(filename);
	if (!in_f.is_open()) {
		std::cerr << std::endl << "ERROR: Unable to open input file." << std::endl;
		return false;
	}

	std::string line;
	key_val_pair pair;
	size_t delim_location;
	size_t line_num = 0;
	while (getNextLine(in_f, line)) {
		++line_num;
		if (line.empty() || line[0] == comment_char) continue;
		if ((delim_location = line.find(delim_char)) != line.rfind(delim_char)) {
			std::cerr << std::endl << "WARN: ignoring line " << line_num << " in input file due to multiple delimiters." << std::endl;
		}
		else if (delim_location == std::string::npos) {
			std::cerr << std::endl << "WARN: ignoring line " << line_num << " in input file due to no delimeter." << std::endl;
		}
		else {
			pair.first = trim(line.substr(0, delim_location));
			pair.first = toupper(pair.first);

			pair.second = trim(line.substr(delim_location + 1));

			if (pair.first.empty()) {
				std::cerr << std::endl << "WARN: ignoring line " << line_num << " in input file due to empty parameter name." << std::endl;
				continue;
			}

			if (pair.second.empty()) {
				std::cerr << std::endl << "WARN: ignoring line " << line_num << " in input file due to empty value." << std::endl;
				continue;
			}

			if (params.find(pair.first) != params.end()) {
				std::cerr << std::endl << "WARN: ignoring line " << line_num << " in input file due to duplicate parameter." << std::endl;
				continue;
			}
			params.insert(pair);

		}
	}

	return true;
}
