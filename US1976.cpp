/*******************************************************************
*   US1976.cpp
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
// US 1976 Atmospheric Density model table
// NOTE: this means the full model has not been implemented here.
// Instead, the predicted density has been tabulated for every
// 1 km of altitude from 0 to 1000 km, and a function is
// provided to linearly interpolate for any requested altitude
// in between.
//

#include "US1976.h"

double get1976Density(const double alt) {
	if (alt <= 0.0) return tbl[0];
	if (alt >= 1e6) return 0.0;
	double alt_below = M_PER_KM*floor(KM_PER_M * alt);
	double alt_above = alt_below + M_PER_KM;
	double rho_below = tbl[static_cast<int>(alt_below) / INT_M_PER_KM];
	double rho_above = tbl[static_cast<int>(alt_above) / INT_M_PER_KM];

	// linearly interpolate between two closest table entries
	// (which are spaced 1 km apart and start at 0 km,
	// i.e. tbl[0] is 0 km, tbl[1] is 1 km, etc.)
	return rho_below + (alt - alt_below)*KM_PER_M*(rho_above - rho_below);
}