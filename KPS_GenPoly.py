#!/usr/bin/env python3.4

####################################################################
#   KPS_GenPoly.py
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
#   Generates a polygon configuration file for input to the KPS
#   propagation framework. The idea is that users can modify this script,
#   or copies of it, to be arbitrarily complex, allowing easy modification
#   of a parameter such as a solar panel deployment angle without having to
#   manually recompute numerical values for all vertices.
#
#   For example, this script currently writes a 6U "dart" configuration
#   satellite, and rather than having to manually compute vertices based
#   on the panel length and angle, simply modify the 'H' and 'PA' variables
#   and run this script to generate the updated polygon file.
#   
#   Verify correct polygon generation with KPS_PlotPoly.
#

from math import sin, cos, pi
from sys import argv

num_vtx = 4

#  folding panel length [m]
H = 0.3

#  panel angle
PA = 130.0 / 180.0 * pi

# header string to go at the top of the output polygon file
header = '# 6U Dart\n# Panel Length: ' + str(H) + '\n# Panel Angle: ' + str(180/pi*PA)

#  The actual body frame satellite as a series of polygons.
#  This is a 6U 'dart' configuration facing in the +z with all 4 panels
#  bilaterally symmetrically deployed to the panel angle specified above.
poly = [\
[0.1, -0.05, -0.15],
[0.1, 0.05, -0.15],
[0.1, 0.05, 0.15],
[0.1, -0.05, 0.15],

[-0.1, -0.05, -0.15],
[-0.1, 0.05, -0.15],
[-0.1, 0.05, 0.15],
[-0.1, -0.05, 0.15],

[-0.1, 0.05, -0.15],
[0.1, 0.05, -0.15],
[0.1, 0.05, 0.15],
[-0.1, 0.05, 0.15],

[-0.1, -0.05, -0.15],
[0.1, -0.05, -0.15],
[0.1, -0.05, 0.15],
[-0.1, -0.05, 0.15],

[-0.1, -0.05, 0.15],
[0.1, -0.05, 0.15],
[0.1, 0.05, 0.15],
[-0.1, 0.05, 0.15],

[-0.1, -0.05, -0.15],
[0.1, -0.05, -0.15],
[0.1, 0.05, -0.15],
[-0.1, 0.05, -0.15],

[0.1, 0.05, -0.15],
[0.1, -0.05, -0.15],
[0.1 + H * sin(PA), -0.05, -0.15 + H * cos(PA)],
[0.1 + H * sin(PA), 0.05, -0.15 + H * cos(PA)],

[-0.1, 0.05, -0.15],
[-0.1, -0.05, -0.15],
[-0.1 - H * sin(PA), -0.05, -0.15 + H * cos(PA)],
[-0.1 - H * sin(PA), 0.05, -0.15 + H * cos(PA)],

[0.1, 0.05, -0.15],
[-0.1, 0.05, -0.15],
[-0.1, 0.05 + H * sin(PA), -0.15 + H * cos(PA)],
[0.1, 0.05 + H * sin(PA), -0.15 + H * cos(PA)],

[0.1, -0.05, -0.15],
[-0.1, -0.05, -0.15],
[-0.1, -0.05 - H * sin(PA), -0.15 + H * cos(PA)],
[0.1, -0.05 - H * sin(PA), -0.15 + H * cos(PA)]]


if len(argv) != 2:
    print('\nUsage: KPS_GenPoly.py <outfile>\n')
else:
    # write polygons to file
    with open(argv[1], 'w') as f:
        f.writelines(tuple(header) + tuple('\n') + tuple(('' if i % num_vtx else '\n') + ', '.join(str(x) for x in poly[i]) + '\n' for i in range(len(poly))))
    print('\nSuccessfully wrote ' + argv[1] + '.\n')
