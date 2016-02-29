#!/usr/bin/env python3.4

####################################################################
#   KPS_GenOrbit.py
#   KPS
#	
#	Author: Kareem Omar
#	kareem.omar@uah.edu
#	https://github.com/komrad36
#
#	Last updated Feb 27, 2016
#   This application is entirely my own work.
####################################################################
#
#   Generates a satellite initial position and velocity
#   in the format required by KPS such that the satellite
#   is in an x by y kilometer orbit at i degrees inclination, where
#   the user specifies x, y, and i. The satellite begins at the 'x'
#   portion of the orbit.
#

from sys import argv, exit
from math import sqrt, sin, cos, pi

if len(argv) != 4:
    print('\nUsage: KPS_GenOrbit.py <x> <y> <inclination>\n\n' \
        'This will generate an <x> by <y> orbit, where x and y are\n' \
        'in km and inclination is in degrees.\n')
    exit(-1)

# Earth standard gravitational parameter
GM = 398600441800000.0

# Earth radius [m]
R  = 6371000.0

M_PER_KM = 1000.0

min_alt = 200
max_alt = 1000

def deg_to_rad(x):
    return x*pi/180.0

x, y, inc = (float(arg) for arg in argv[1:])

if x < min_alt or y < min_alt:
    print('ERROR: Altitudes must be above', min_alt, 'km. Aborting.\n')
    exit(-1)

if x > max_alt or y > max_alt:
    print('WARNING: KPS is intended for altitudes below', max_alt, \
       'km.\nProceed at your own risk.\n')

if inc < -90.0 or inc > 90.0:
    print('ERROR: inclination must be -90 to 90', chr(176), '. Aborting.\n', sep='')
    exit(-1)

x_r = R + M_PER_KM*x
y_r = R + M_PER_KM*y

init_pos = [x_r, 0.0, 0.0]

# semi-major axis
a = 0.5*(x_r + y_r)

# vis-viva
init_vel_mag = sqrt(GM*(2/x_r-1/a))
init_vel = [0.0, init_vel_mag*cos(deg_to_rad(inc)), init_vel_mag*sin(deg_to_rad(inc))]

print('\nGenerating stats for a {0:.14g}x{1:.14g} km orbit at {2:.14g}{3}...\n'.format(x, y, inc, chr(176)))
print('SAT_INIT_POS = {0:.14g}, {1:.14g}, {2:.14g}'.format(*init_pos))
print('SAT_INIT_V = {0:.14g}, {1:.14g}, {2:.14g}\n'.format(*init_vel))
