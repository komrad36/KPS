/*******************************************************************
*   PolyInput.h
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

#ifndef POLY_INPUT_H
#define POLY_INPUT_H

#include <algorithm>
#include <fstream>
#include <iostream>

#include "FileInput.h"
#include "Params.h"
#include "str_util.h"

const int DIMS = 3;
const int MAX_SUPPORTED_POLY = 512;

bool parsePoly(std::vector<vec3>& poly, int& num_poly, const std::string& filename, const char comment_char = '#');

#endif
