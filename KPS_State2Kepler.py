#!/usr/bin/env python3.4

####################################################################
#   KPS_State2Kepler.py
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
#   Converts state vectors to Keplerian orbital elements.
#
#   A state vector may be directly supplied, OR a single number
#   can be supplied - a time in seconds. KPS_State2Kepler will
#   search through the most recent KPS run and retrieve the
#   first state vector AFTER that time for conversion.
#

from sys import argv
from math import sin, cos, tan, acos, atan2, ceil, pi, sqrt

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

if len(argv) != 2 and len(argv) != 7:
        print('\nUsage: KPS_State2Kepler.py <time in seconds>\n'
          'to read data for <time> from an existing KPS run\n'
          '  -- OR --\n'
          'KPS_State2Kepler.py <rx> <ry> <rz> <vx> <vy> <vz>\n'
          'where r and v are position and velocity in inertial frame, in m and m/s\n')
else:
    # if in time retrieval mode
    if len(argv) == 2:
        # grab data from files
        import numpy as np
        from sys import exit
        try:
            f_t = open('t.bin', 'rb')
        except:
            print('ERROR: failed to open file. Aborting.')
            exit(-1)

        # discard last entry in case it is corrupt
        # (if it's being written realtime)
        t = np.fromfile(f_t)[:-1]

        try:
            f_r = open('r.bin', 'rb')
        except:
            print('ERROR: failed to open file. Aborting.')
            exit(-1)
        r = np.fromfile(f_r)
        r = r[:3*ceil(r.size//3) - 2]

        try:
            f_v = open('v.bin', 'rb')
        except:
            print('ERROR: failed to open file. Aborting.')
            exit(-1)
        v = np.fromfile(f_v)
        v = v[:3*ceil(v.size//3) - 2]

        f_t.close()
        f_r.close()
        f_v.close()

        t_val = float(argv[1])

        if t_val < t[0]:
            print('ERROR: requested time not in range of files! Aborting.')
            exit(-1)

        for i in range(len(t)):
            if t[i] >= t_val:
                # found the first time step after requested time
                r = r[3*i:3*i+3].tolist()
                v = v[3*i:3*i+3].tolist()
                print('\nTime: {0:.14g} s'.format(t[i]))
                del t
                break
        else:
            print('ERROR: requested time not in range of files! Aborting.')
            exit(-1)
    else:
        r = [float(val) for val in argv[1:4]]
        v = [float(val) for val in argv[4:]]

    # Earth gravitational parameter
    GM = 398600441800000.0

    v_mag = norm(v)
    r_mag = norm(r)

    # orbital momentum
    h = cross(r, v)

    # eccentricity vector
    e_vec = [cross(v, h)[i]/GM - r[i]/r_mag for i in range(3)]
    
    # ascending node vector
    n = [-h[1], h[0], 0.0]
    n_mag = norm(n)

    # true anomaly
    nu = acos(dot(e_vec, r)/(norm(e_vec)*r_mag))
    if dot(r, v) < 0.0:
        nu = 2*pi - nu

    # inclination
    i = acos(h[2]/norm(h))

    # scalar eccentricity
    e = norm(e_vec)

    # eccentric anomaly
    E = 2*atan2(tan(0.5*nu), sqrt((1+e)/(1-e)))

    # RAAN
    Omega = acos(n[0]/n_mag)
    if n[1] < 0.0:
        Omega = 2*pi - Omega

    # argument of periapsis
    w = acos(dot(n, e_vec)/(n_mag*e))
    if e_vec[2] < 0.0:
        w = 2*pi - w

    # mean anomaly
    M = E - e*sin(E)

    # semi-major axis
    a = 1.0/(2.0/r_mag-v_mag*v_mag/GM)

    nu = rad_to_deg(nu)
    i = rad_to_deg(i)
    E = rad_to_deg(E)
    Omega = rad_to_deg(Omega)
    w = rad_to_deg(w)
    M = rad_to_deg(M)

    print('\nPosition Vector: {0:.14g}, {1:.14g}, {2:.14g} m'.format(r[0], r[1], r[2]))
    print('Velocity Vector: {0:.14g}, {1:.14g}, {2:.14g} m/s'.format(v[0], v[1], v[2]))

    print('True Anomaly: {0:.14g}{1}'.format(nu, chr(176)))
    print('Mean Anomaly: {0:.14g}{1}'.format(M, chr(176)))
    print('Eccentric Anomaly: {0:.14g}{1}'.format(E, chr(176)))

    print('Semi-major Axis: {0:.14g} m'.format(a))
    print('Eccentricity: {0:.14g}'.format(e))
    print('Argument of Periapsis: {0:.14g}{1}'.format(w, chr(176)))
    print('Longitude of Ascending Node: {0:.14g}{1}'.format(Omega, chr(176)))
    print('Inclination: {0:.14g}{1}\n'.format(i, chr(176)))
