/*******************************************************************
*   ParamInput.h
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
// Parses a KPS parameters file into a map of parameter names and values.
//

#ifndef PARAM_INPUT_H
#define PARAM_INPUT_H

#include <fstream>
#include <iostream>

#include "FileInput.h"
#include "Params.h"
#include "str_util.h"

// parses 'filename' into 'params', a map of parameter names and values, ignoring lines that begin with
// 'comment_char' and empty lines and using 'delim_char' to separate param names and values
// DEFAULT PARAMETER: delim_char = '='
// DEFAULT PARAMETER: comment_char = '#'
bool parseParams(key_val_map& params, const std::string& filename, const char delim_char = '=', const char comment_char = '#');

#endif
