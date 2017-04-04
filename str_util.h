/*******************************************************************
*   str_util.h
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
// Includes and inline convenience functions for manipulating
// strings that are not already available in the standard library
//

#pragma once

#include <algorithm>
#include <cctype>
#include <functional>
#include <string>

// trim whitespace from head and tail of string
// empty strings are handled safely
inline std::string trim(const std::string &s) {
	auto front = std::find_if_not(s.begin(), s.end(), ::isspace);
	auto back = std::find_if_not(s.rbegin(), s.rend(), ::isspace).base();
	return (back <= front ? std::string() : std::string(front, back));
}

// tokenize a string into a vector of substrings (without the delimiter)
// DEFAULT PARAMETER delim = ','
inline void tokenize(std::vector<std::string>& vals, const std::string& s, const char delim = ',') {
	vals.clear();
	size_t start = 0;
	for (size_t i = 0;;) {
		if (i >= s.length()) break;
		if (s[i++] == delim) {
			vals.push_back(trim(s.substr(start, i - start - 1)).c_str());
			start = i;
		}
	}
	vals.push_back(trim(s.substr(start).c_str()));
}

// transform an entire string to lowercase, ignoring non-alpha characters
inline std::string tolower(std::string s) {
	std::transform(s.begin(), s.end(), s.begin(), std::ptr_fun<int, int>(std::tolower));
	return s;
}

// transform an entire string to uppercase, ignoring non-alpha characters
inline std::string toupper(std::string s) {
	std::transform(s.begin(), s.end(), s.begin(), std::ptr_fun<int, int>(std::toupper));
	return s;
}


