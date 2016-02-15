/*******************************************************************
*   PolyInput.cpp
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
// Parses a KPS polygon file into a vector of vertices.
//

#include "PolyInput.h"

// parse 'filename' into polygons stored in 'poly' as a simple, flat vector of vertices
// DEFAULT PARAMETER: comment_char = '#'
bool parsePoly(std::vector<vec3>& poly, int& num_poly, const std::string& filename, const char comment_char) {
	std::ifstream in_f(filename);
	if (!in_f.is_open()) {
		std::cerr << std::endl << "ERROR: Unable to open input file." << std::endl;
		return false;
	}

	std::string line;
	size_t line_num = 0;
	std::vector<std::string> vals;
	poly.clear();
	while (getNextLine(in_f, line)) {
		++line_num;
		if (line.empty() || line[0] == comment_char) continue;

		tokenize(vals, line);
		if (vals.size() != DIMS) {
			std::cerr << std::endl << "ERROR: line " << line_num << " of polygon file has invalid vertex. Aborting." << std::endl;
			return false;
		}

		vec3 vtx;
		for (int i = 0; i < DIMS; ++i) {
			vtx[i] = atof(vals[i].c_str());
		}
		poly.push_back(vtx);
	}

	if (poly.size() % NUM_VTX) {
		std::cerr << std::endl << "ERROR: " << poly.size() << " vertices read. KPS is currently compiled to expect" << std::endl
			<< NUM_VTX << " vertices per polygon, so the polygon file has incomplete polygons. Aborting." << std::endl;
		return false;
	}

	num_poly = static_cast<int>(poly.size()) / NUM_VTX;

	if (num_poly > MAX_SUPPORTED_POLY) {
		std::cerr << std::endl << "ERROR: " << num_poly << " polygons loaded, but KPS is currently compiled to support a" << std::endl
			<< "maximum of " << MAX_SUPPORTED_POLY << ". Aborting." << std::endl;
		return false;
	}

	// all points in a polygon must be coplanar, but due to truncation error,
	// 16 times machine epsilon is a reasonable tolerance to check against
	const double tol = 16.0 * std::numeric_limits<double>::epsilon();
	for (int i = 0; i < num_poly; ++i) {

		vec3 normal = glm::cross(poly[i*NUM_VTX + 1] - poly[i*NUM_VTX], poly[i*NUM_VTX + 2] - poly[i*NUM_VTX]);
		for (int j = i*NUM_VTX + 3; j < (i + 1)*NUM_VTX; ++j) {
			if (fabs(glm::dot(poly[j] - poly[i*NUM_VTX], normal)) > tol) {
				std::cerr << std::endl << "ERROR: non-planar polygon found! Aborting." << std::endl;
				return false;
			}
		}

		for (int j = i*NUM_VTX; j < (i+1)*NUM_VTX; ++j) {
			for (int k = i*NUM_VTX; k < (i+1)*NUM_VTX; ++k) {
				if (j != k && poly[j] == poly[k]) {
					std::cerr << std::endl << "ERROR: duplicate vertex found within a single polygon. Aborting." << std::endl;
					return false;
				}
			}
		}

	}

	poly.shrink_to_fit();
	return true;
}
