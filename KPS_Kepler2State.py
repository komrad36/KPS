#!/usr/bin/env python3.4

####################################################################
#   KPS_Kepler2State.py
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
#   Converts Keplerian orbital elements to state vector.
#         Inputs:
#         a      -  Semi-major Axis [m]
#         e      -  Eccentricity
#         i      -  Inclination [degrees]
#         Omega  -  Longitude of Ascending Node [degrees]
#         w      -  Argument of Periapsis [degrees]
#         M      -  Mean Anomaly [degrees]
#

from sys import argv, exit
from math import sin, cos, tan, acos, atan2, ceil, pi, sqrt

MAX_ITER = 100000

def norm(vec):
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2])

def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]]

def dot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def rad_to_deg(x):
    return x * 180.0 / pi

def deg_to_rad(x):
    return x * pi / 180.0

if len(argv) != 7:
        print('\nUsage: KPS_Kepler2State.py <a> <e> <i> <Omega> <w> <M>\n'
          ' Inputs:\n'
          '   a      -  Semi-major Axis [m]\n'
          '   e      -  Eccentricity\n'
          '   i      -  Inclination [degrees]\n'
          '   Omega  -  Longitude of Ascending Node [degrees]\n'
          '   w      -  Argument of Periapsis [degrees]\n'
          '   M      -  Mean Anomaly [degrees]\n')
else:
    [a, e, i, Omega, w, M] = (float(val) for val in argv[1:])

    # Earth standard gravitational parameter
    GM = 398600441800000.0

    w = deg_to_rad(w)
    Omega = deg_to_rad(Omega)
    i = deg_to_rad(i)
    M = deg_to_rad(M)

    # Newton's method solve of Kepler's equation for eccentric anomaly
    E_last, E, iter = None, M, 0
    while(E_last != E):
        if iter > MAX_ITER:
            print('ERROR: solution does not converge! Check inputs and try again.')
            exit(-1)
        E_last = E
        E = E_last - (E_last-e*sin(E_last)-M)/(1-e*cos(E_last))
        iter += 1

    # true anomaly
    nu = 2*atan2(sqrt(1+e)*sin(0.5*E), sqrt(1-e)*cos(0.5*E))

    # distance to Earth center
    r_c = a*(1-e*cos(E))

    # orbital frame
    r_o = [r_c*cos(nu), r_c*sin(nu), 0.0]
    v_o = [-sqrt(GM*a)/r_c*sin(E), sqrt(GM*a)/r_c*sqrt(1-e*e)*cos(E), 0.0]

    # to ECI frame
    r = [r_o[0]*(cos(w)*cos(Omega)-sin(w)*cos(i)*sin(Omega)) - r_o[1]*(sin(w)*cos(Omega)+cos(w)*cos(i)*sin(Omega)),
        r_o[0]*(cos(w)*sin(Omega)+sin(w)*cos(i)*cos(Omega)) + r_o[1]*(cos(w)*cos(i)*cos(Omega)-sin(w)*sin(Omega)),
        r_o[0]*sin(w)*sin(i) + r_o[1]*cos(w)*sin(i)]

    v = [v_o[0]*(cos(w)*cos(Omega)-sin(w)*cos(i)*sin(Omega)) - v_o[1]*(sin(w)*cos(Omega)+cos(w)*cos(i)*sin(Omega)),
        v_o[0]*(cos(w)*sin(Omega)+sin(w)*cos(i)*cos(Omega)) + v_o[1]*(cos(w)*cos(i)*cos(Omega)-sin(w)*sin(Omega)),
        v_o[0]*sin(w)*sin(i) + v_o[1]*cos(w)*sin(i)]

    w = rad_to_deg(w)
    Omega = rad_to_deg(Omega)
    i = rad_to_deg(i)
    nu = rad_to_deg(nu)
    M = rad_to_deg(M)
    E = rad_to_deg(E)

    print('\nSemi-major Axis: {0:.14g} m'.format(a))
    print('Eccentricity: {0:.14g}'.format(e))
    print('Argument of Periapsis: {0:.14g}{1}'.format(w, chr(176)))
    print('Longitude of Ascending Node: {0:.14g}{1}'.format(Omega, chr(176)))
    print('Inclination: {0:.14g}{1}'.format(i, chr(176)))

    print('True Anomaly: {0:.14g}{1}'.format(nu, chr(176)))
    print('Mean Anomaly: {0:.14g}{1}'.format(M, chr(176)))
    print('Eccentric Anomaly: {0:.14g}{1}'.format(E, chr(176)))

    print('Position Vector: {0:.14g}, {1:.14g}, {2:.14g} m'.format(r[0], r[1], r[2]))
    print('Velocity Vector: {0:.14g}, {1:.14g}, {2:.14g} m/s\n'.format(v[0], v[1], v[2]))
