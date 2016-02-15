/*******************************************************************
*   FileInput.h
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
// Platform-indepdendent file reading (i.e. independent of newline character)
//

#ifndef FILEIO_H
#define FILEIO_H

#include <fstream>
#include <string>

#define DO_NOT_SKIP_WHITESPACE	(true)

bool getNextLine(std::ifstream& in_f, std::string& line);

#endif