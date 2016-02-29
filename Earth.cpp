/*******************************************************************
*   Earth.cpp
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
// Earth (orbited body) model, including wrappers for gravitational modeling
//

#include "Earth.h"

// standard gravitational parameter
const double Earth::OMEGA = 7.29211509e-5;

// rate of Earth's rotation about its axis [rad/sec]
const double Earth::GM = 398600441800000.0;

// Earth's radius [m]
const double Earth::R = 6371000.0;

// altitude at which deorbit should be assumed [m]
const double Earth::DEORBIT_ALT = 150000.0;
