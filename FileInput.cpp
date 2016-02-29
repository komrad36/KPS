/*******************************************************************
*   FileInput.cpp
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
// Platform-indepdendent file reading (i.e. independent of newline character)
//

#include "FileInput.h"

// hand-rolled getLine compatible with
// Windows and Linux newlines, as well
// as files without a final trailing newline
std::istream& safeGetline(std::istream& is, std::string& t) {
	t.clear();

	// streambuf for speed
	// requires stream sentry to guard buffer
	std::istream::sentry se(is, DO_NOT_SKIP_WHITESPACE);
	std::streambuf* sb = is.rdbuf();

	short c;

	for (;;) {
		c = sb->sbumpc();
		switch (c) {
		case '\n':
			return is;
		case '\r':
			if (sb->sgetc() == '\n')
				sb->sbumpc();
			return is;
		case EOF:
			// handle the case when the last line has no line ending
			if (t.empty())
				is.setstate(std::ios::eofbit);
			return is;
		default:
			t += static_cast<char>(c);
		}
	}
}

// wrapper to simply retrive clean line from ifstream
// returns false when no more valid lines
bool getNextLine(std::ifstream& in_f, std::string& line) {
	bool cont = true;

	while (cont) {
		safeGetline(in_f, line);

		if (!in_f.eof()) {
			cont = false;
		}
		else {
			line.clear();
			return false;
		}

	}

	return true;
}