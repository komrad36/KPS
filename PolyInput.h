/*******************************************************************
*   PolyInput.h
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
// Parses a KPS polygon file into a vector of vertices.
//

#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>

#include "FileInput.h"
#include "Params.h"
#include "str_util.h"

constexpr int DIMS = 3;
constexpr int MAX_SUPPORTED_POLY = 512;

bool parsePoly(std::vector<vec3>& poly, int& num_poly, const std::string& filename, const char comment_char = '#');


