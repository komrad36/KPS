#!/usr/bin/env python3.4

####################################################################
#   KPS_PlotPoly.pyw
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
#   Plots polygons of generated polygon file in 3-D to
#   examine and verify correctness.
#

font_size = 21
maximize_plot = True
wireframe = False
poly_color = 'blue'
num_vtx = 4
face_alpha = 0.75

# maximize plots if desired, on any backend
def maximizePlot():
    try:
        mng = plt.get_current_fig_manager()
        backend = plt.get_backend()
        if backend == 'TkAgg':
            try:
                mng.window.state('zoomed')
            except:
                mng.resize(*mng.window.maxsize())
        elif backend == 'wxAgg':
            mng.frame.Maximize(True)
        elif backend[:2].upper() == 'QT':
            mng.window.showMaximized()
        else:
            return False
        return True
    except:
        return False

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D
import sys

if len(sys.argv) != 2:
    print('\nUsage: KPS_PlotPoly.py <polygon_file>\n')
else:
    name = sys.argv[1]
    short_name = name[max(name.rfind('/'), name.rfind('\\')) + 1:]
    try:
        f = open(name)
    except:
        print('Failed to open ' + name + '. Aborting.')

    plt.rcParams['font.size'] = font_size

    fig = plt.figure('KPS - ' + short_name)
    ax = fig.gca(projection='3d', aspect='equal')

    # load polygons
    vtx = [[float(x) for x in line.rstrip('\n').split(',')] for line in f if len(line) > 4 and line[0] != '#']
    x, y, z = list(map(list, zip(*vtx)))
    poly = [vtx[i:i + num_vtx] for i in range(0, len(vtx), num_vtx)]

    ax.set_xlabel(r'$x_{body}  [m]$')
    ax.set_ylabel(r'$y_{body}  [m]$')
    ax.set_zlabel(r'$z_{body}  [m]$')

    x_min = min(x)
    x_max = max(x)
    x_center = 0.5 * (x_max + x_min)

    y_min = min(y)
    y_max = max(y)
    y_center = 0.5 * (y_max + y_min)

    z_min = min(z)
    z_max = max(z)
    z_center = 0.5 * (z_max + z_min)

    total_min = min([x_min, y_min, z_min])
    total_max = max([x_max, y_max, z_max])
    half_span = 0.5 * (total_max - total_min)

    ax.set_xlim3d(x_center - half_span, x_center + half_span)
    ax.set_ylim3d(y_center - half_span, y_center + half_span)
    ax.set_zlim3d(z_center - half_span, z_center + half_span)

    if wireframe:
        for i in range(0, len(z), num_vtx):
            ax.plot(x[i:i + num_vtx] + [x[i]], y[i:i + num_vtx] + [y[i]], z[i:i + num_vtx] + [z[i]], color=poly_color)
    else:
        ax.add_collection3d(Poly3DCollection(poly, alpha=face_alpha, edgecolor='k', color=poly_color))

    ax.plot([x_min, x_max], [0, 0], [0, 0], color='k', alpha=0.35)
    ax.plot([0, 0], [y_min, y_max], [0, 0], color='k', alpha=0.35)
    ax.plot([0, 0], [0, 0], [z_min, z_max], color='k', alpha=0.35)

    fig.tight_layout()

    if maximize_plot:
        maximizePlot()
    fig.canvas.draw()
    plt.show()
