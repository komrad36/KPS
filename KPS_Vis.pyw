#!/usr/bin/env python3.4

####################################################################
#   KPS_Vis.py
#   KPS
#	
#	Author: Kareem Omar
#	kareem.omar@uah.edu
#	https://github.com/komrad36
#
#	Last updated Mar 20, 2016
#   This application is entirely my own work.
####################################################################
#
#   Realtime visualization of any combination of the eight parameters
#   selectable below. Just comment out the ones you don't want.
#
#   Visualization begins in realtime mode, constantly checking
#   the outfiles for updates. When no changes occur in the files
#   for some time, KPS_Vis switches to FINAL mode, adjusting the axes
#   and allowing user interaction.
#

# do not modify this line; configure plotting options below
t, R, V, Q, Q_ORB, ALT, B_STAR, W, V_B, E, ORIENTATION, SEMI_MAJOR, ECC, INC, RAAN, PERIAPSIS, MEAN_ANOM, TRUE_ANOM, ECC_ANOM = range(19)

#####################################################################
##################  USER CONFIGURABLE PARAMETERS  ###################
#####################################################################

font_size = 14
line_width = 1.5
maximize_plot = True

# for live orientation plotter
wireframe = False
poly_color = 'green'
num_vtx = 4
poly_file = 'poly.kps'
face_alpha = 0.8
axis_alpha = 0.5

# simply comment
# out the ones you don't want
plots = [
#R,
ORIENTATION,
#V,
#Q,
#Q_ORB,
#ALT,
#B_STAR,
#W,
#V_B,
E,
#SEMI_MAJOR,
#ECC,
#INC,
#RAAN,
#PERIAPSIS,
#MEAN_ANOM,
#TRUE_ANOM,
#ECC_ANOM,
]

#####################################################################
#####################################################################
#####################################################################




from math import pi, ceil, sqrt, acos, atan2, sin, tan
import matplotlib.pyplot as plt
import numpy as np
from tkinter import TclError

if any(x == R or x == ORIENTATION for x in plots):
    from mpl_toolkits.mplot3d import Axes3D

if any(x == ORIENTATION for x in plots):
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]]

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
        elif backend == 'QT4Agg':
            mng.window.showMaximized()
        else:
            return False
        return True
    except:
        return False

# number of elements for each input
elems = [1,3,3,4,4,1,1,3,3,1,0,1,1,1,1,1,1,1,1]

# names of data files
names = ['t', 'r', 'v', 'q', 'q_orb', 'alt', 'b_star', 'w', 'v_body', 'p_e']
titles = ['', 'Position', 'Velocity', 'Quaternion to ECI Frame', 'Quaternion to Orbital Frame', 'Altitude', 'Starred Ballistic Coefficient', 'Angular Velocity', 'Velocity in Body Frame', 'Pointing Error', 'Live Orientation', 'Semi-Major Axis', 'Eccentricity', 'Inclination', 'Longitude of Ascending Node', 'Argument of Periapsis', 'Mean Anomaly', 'True Anomaly', 'Eccentric Anomaly']
ylabels = [r'', r'$y_{ECI}  [m]$', '[m/s]', 'Component', 'Component', '[m]', r'[$m^{-1}$]', r'$[s^{-1}]$', '[m/s]', '['+chr(176)+']', r'$y_{orbital}  [m]$', '[m]', [], '['+chr(176)+']', '['+chr(176)+']', '['+chr(176)+']', '['+chr(176)+']', '['+chr(176)+']', '['+chr(176)+']']
xlabels = [r'$x_{ECI}  [m]$', r'$x_{orbital}  [m]$']
zlabels = [r'$z_{ECI}  [m]$', r'$z_{orbital}  [m]$']
num_plots = len(plots)
compute_keplerian = False

# data files required by the user's choice of plots
reqd = [t]
for j in range(num_plots):
    i = plots[j]
    if i < ORIENTATION and not any(x == i for x in reqd):
        reqd.append(i)
    elif i == ORIENTATION and not any(x == Q_ORB for x in reqd):
        reqd.append(Q_ORB)
        try:
            f_orb = open(poly_file)
        except:
            print('ERROR: failed to open polygon file. Aborting.')
            sys.exit(-1)
        vtx = np.array([[float(x) for x in line.rstrip('\n').split(', ')] for line in f_orb if len(line) > 4 and line[0] != '#'])
        num_poly = len(vtx) // num_vtx
        f_orb.close()
    else:
        if not any(x == R for x in reqd):
            reqd.append(R)
        if not any(x == V for x in reqd):
            reqd.append(V)
        compute_keplerian = True

num_data = len(reqd)
h_f = [[] for _ in range(19)]
data = [[] for _ in range(19)]

for j in range(num_data):
    i = reqd[j]
    try:
        h_f[i] = open(names[i] + '.bin', 'rb')
    except:
        print('ERROR: failed to open file. Aborting.')
        sys.exit(-1)

plt.rcParams.update({'font.size': font_size})
h_fig = plt.figure('KPS - REALTIME')

if maximize_plot:
    maximizePlot()

plt.show(False)

f_blank = [[],
           lambda: [plt.plot([], [], [], color='g', linewidth=line_width)[0],
                    plt.plot([0], [0], [0], 'bo', markersize=40, alpha=0.90)[0]],
           lambda: [plt.plot([], [], color='r', label=r'$v_x$', linewidth=line_width)[0],
                    plt.plot([], [], color='g', label=r'$v_y$', linewidth=line_width)[0],
                    plt.plot([], [], color='b', label=r'$v_z$', linewidth=line_width)[0]],
           lambda: [plt.plot([], [], color='r', label=r'$q_0$', linewidth=line_width)[0],
                    plt.plot([], [], color='g', label=r'$q_1$', linewidth=line_width)[0],
                    plt.plot([], [], color='b', label=r'$q_2$', linewidth=line_width)[0],
                    plt.plot([], [], color='k', label=r'$q_3$', linewidth=line_width)[0]],
           lambda: [plt.plot([], [], color='r', label=r'$q_0$', linewidth=line_width)[0],
                    plt.plot([], [], color='g', label=r'$q_1$', linewidth=line_width)[0],
                    plt.plot([], [], color='b', label=r'$q_2$', linewidth=line_width)[0],
                    plt.plot([], [], color='k', label=r'$q_3$', linewidth=line_width)[0]],
           lambda: [plt.plot([], [], color='g', label=r'alt', linewidth=line_width)[0]],
           lambda: [plt.plot([], [], color='b', label=r'$B^*$', linewidth=line_width)[0]],
           lambda: [plt.plot([], [], color='r', label=r'$\omega_x$', linewidth=line_width)[0],
                    plt.plot([], [], color='g', label=r'$\omega_y$', linewidth=line_width)[0],
                    plt.plot([], [], color='b', label=r'$\omega_z$', linewidth=line_width)[0]],
           lambda: [plt.plot([], [], color='r', label=r'$v_x$', linewidth=line_width)[0],
                    plt.plot([], [], color='g', label=r'$v_y$', linewidth=line_width)[0],
                    plt.plot([], [], color='b', label=r'$v_z$', linewidth=line_width)[0]],
           lambda: [plt.plot([], [], color='r', label=r'$\epsilon_p$', linewidth=line_width)[0]],
           lambda: [],
           lambda: [plt.plot([], [], color='g', label=r'a', linewidth=line_width)[0]],
           lambda: [plt.plot([], [], color='b', label=r'e', linewidth=line_width)[0]],
           lambda: [plt.plot([], [], color='r', label=r'i', linewidth=line_width)[0]],
           lambda: [plt.plot([], [], color='g', label=r'$\Omega$', linewidth=line_width)[0]],
           lambda: [plt.plot([], [], color='b', label=r'$\omega$', linewidth=line_width)[0]],
           lambda: [plt.plot([], [], color='r', label=r'M', linewidth=line_width)[0]],
           lambda: [plt.plot([], [], color='g', label=r'$\nu$', linewidth=line_width)[0]],
           lambda: [plt.plot([], [], color='b', label=r'E', linewidth=line_width)[0]]]

# init plots
h_plots = [[] for _ in range(19)]
h_lines = [[] for _ in range(19)]
for j in range(num_plots):
    i = plots[j]
    if i == R or i == ORIENTATION:
        h_plots[i] = h_fig.add_subplot(num_plots, 1, j + 1, aspect='equal', projection='3d')
    else:
        h_plots[i] = h_fig.add_subplot(num_plots, 1, j + 1)
 
    if i == ORIENTATION:
        h_lines[i].append([])
        if wireframe:
            for k in range(num_poly):
                h_lines[i][0].append(h_plots[i].plot([], [], [], color=poly_color, linewidth=line_width)[0])
        else:
            p_collection = Poly3DCollection([], alpha=face_alpha, edgecolor='k', color=poly_color)
            h_plots[i].add_collection3d(p_collection)

        h_lines[i].append(h_plots[i].plot([], [], [], color='k', alpha=axis_alpha)[0])
        h_lines[i].append(h_plots[i].plot([], [], [], color='k', alpha=axis_alpha)[0])
        h_lines[i].append(h_plots[i].plot([], [], [], color='k', alpha=axis_alpha)[0])
    else:
        h_lines[i] = f_blank[i]()

    if i != R and i != ORIENTATION:
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.title(titles[i], fontsize=1.2*font_size, fontweight="bold")
    plt.ylabel(ylabels[i], fontsize=font_size)
    
    if i == E:
        plt.ylim([0, 180])
    if i == INC:
        plt.ylim([-90, 90])
    if i >= RAAN:
        plt.ylim([0, 360])
    if i == R or i == ORIENTATION:
        plt.xlabel(xlabels[i!=R])
        h_plots[i].set_zlabel(zlabels[i!=R])

# print xlabel if there are plots other than the 3-D ones, orientation and position 
if (num_plots == 1 and (plots[0] == ORIENTATION or plots[0] == R)) or (num_plots == 2 and any(x == ORIENTATION for x in plots) and any(x == R for x in plots)):
    h_fig.tight_layout()
else:
    plt.xlabel('Time [s]', fontsize=1.4*font_size)

h_fig.canvas.draw()
canvas = [[] for _ in range(19)]

def onresize(event):
    global canvas
    for j in range(num_plots):
        i = plots[j]
        canvas[i] = h_fig.canvas.copy_from_bbox(h_plots[i].bbox)

onresize(None)

remove_text = False
if plots and plots[-1] != R and plots[-1] != ORIENTATION:
    remove_text = True
    h_text = h_plots[plots[-1]].text(0.5, 0.5, 'Loading...', fontsize=50, horizontalalignment='center', verticalalignment='center', transform=h_plots[plots[-1]].transAxes)
h_fig.canvas.draw()

ticks_without_change = 0
rad_to_deg = 180.0/pi

continue_please = True

def onclose(event):
    global continue_please
    continue_please = False

# force redraw on resize
cid = h_fig.canvas.mpl_connect('resize_event', onresize)

# graceful close
h_fig.canvas.mpl_connect('close_event', onclose)

try:
    while continue_please:
        old_t_size = len(data[0])
        # update live data
        for j in range(num_data):
            i = reqd[j]
            h_f[i].seek(8*len(data[i]))
            data[i].extend(np.fromfile(h_f[i]))

        # shorten all vectors to length of shortest one
        # minus 1, since with realtime data the very last
        # entry might still have been incomplete and thus
        # corrupt
        min_len = min(len(data[x])//elems[x] for x in reqd) - 1
        for j in range(num_data):
            i = reqd[j]
            data[i] = data[i][:elems[i]*min_len]

        if min_len == old_t_size:
            ticks_without_change += 1
            if ticks_without_change > 20:
                # no change for a while
                # switch to Final mode
                h_fig.canvas.set_window_title('KPS - Final')
                h_fig.canvas.mpl_disconnect(cid)
                for j in range(num_plots):
                    i = plots[j]
                    if i != R and i != ORIENTATION:
                        h_plots[i].set_xlim([0, data[0][-1]])
                        h_plots[i].autoscale_view(True, True, False)
                plt.show()
        else:
            ticks_without_change = 0

            # only compute Keplerian elements if necessary
            if compute_keplerian:
                rng_a = len(data[SEMI_MAJOR])
                delta = min_len - rng_a
                # preallocate new section of memory
                data[SEMI_MAJOR].extend([0]*delta)
                data[TRUE_ANOM].extend([0]*delta)
                data[INC].extend([0]*delta)
                data[ECC].extend([0]*delta)
                data[RAAN].extend([0]*delta)
                data[PERIAPSIS].extend([0]*delta)
                data[MEAN_ANOM].extend([0]*delta)
                data[ECC_ANOM].extend([0]*delta)

                for k in range(rng_a, min_len):
                    r = data[R][3*k:3*k+3]
                    v = data[V][3*k:3*k+3]

                    v_mag = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
                    r_mag = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2])

                    # orbital momentum
                    h = cross(r, v)

                    # eccentricity vector
                    # the magic num is 1/mu, Earth's standard gravitational parameter,
                    # precomputed for speed
                    e_vec = [2.508777951886354381858e-15*cross(v, h)[x] - r[x]/r_mag for x in range(3)]
                    # scalar eccentricity
                    e = sqrt(e_vec[0]*e_vec[0] + e_vec[1]*e_vec[1] + e_vec[2]*e_vec[2])
                
                    # ascending node vector
                    n_1 = -h[1]
                    n_2 = h[0]
                    # n_3 is 0
                    n_mag = sqrt(n_1*n_1 + n_2*n_2)

                    # true anomaly
                    nu = acos((e_vec[0]*r[0]+e_vec[1]*r[1]+e_vec[2]*r[2])/(e*r_mag))
                    if r[0]*v[0] + r[1]*v[1] + r[2]*v[2] < 0.0:
                        nu = 2*pi - nu

                    # inclination
                    inc = acos(h[2]/sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]))
                
                    # eccentric anomaly
                    E = 2*atan2(tan(0.5*nu), sqrt((1+e)/(1-e)))
                    if E < 0.0:
                        E += 2*pi

                    # RAAN
                    Omega = acos(n_1/n_mag)
                    if n_2 < 0.0:
                        Omega = 2*pi - Omega

                    # argument of periapsis
                    w = acos((n_1*e_vec[0] + n_2*e_vec[1])/(n_mag*e))
                    if e_vec[2] < 0.0:
                        w = 2*pi - w

                    # mean anomaly
                    M = E - e*sin(E)
                    if M < 0.0:
                        M += 2*pi

                    # the magic num is 1/mu, Earth's standard gravitational parameter,
                    # precomputed for speed
                    data[SEMI_MAJOR][k] = 1.0/(2.0/r_mag-2.508777951886354381858e-15*v_mag*v_mag)
                    data[TRUE_ANOM][k] = rad_to_deg*nu
                    data[INC][k] = rad_to_deg*inc
                    data[ECC][k] = e
                    data[RAAN][k] = rad_to_deg*Omega
                    data[PERIAPSIS][k] = rad_to_deg*w
                    data[MEAN_ANOM][k] = rad_to_deg*M
                    data[ECC_ANOM][k] = rad_to_deg*E
            
            if remove_text:
                h_text.remove()
                remove_text = False
            need_redraw = False
            # plot updated data
            for j in range(num_plots):
                i = plots[j]
                if i == ORIENTATION:
                    q_vec = data[Q_ORB][-3:]
                    vtx_rot = [v + cross([2*q_vec[0], 2*q_vec[1], 2*q_vec[2]],cross(q_vec, v)+data[Q_ORB][-4]*v) for v in vtx]

                    x, y, z = list(map(list, zip(*vtx_rot)))

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

                    h_plots[ORIENTATION].set_xlim3d(x_center - half_span, x_center + half_span)
                    h_plots[ORIENTATION].set_ylim3d(y_center - half_span, y_center + half_span)
                    h_plots[ORIENTATION].set_zlim3d(z_center - half_span, z_center + half_span)

                    if wireframe:
                        for k in range(0, len(z), num_vtx):
                            h_lines[ORIENTATION][0][k//num_vtx].set_xdata(x[k:k + num_vtx] + [x[k]])
                            h_lines[ORIENTATION][0][k//num_vtx].set_ydata(y[k:k + num_vtx] + [y[k]])
                            h_lines[ORIENTATION][0][k//num_vtx].set_3d_properties(z[k:k + num_vtx] + [z[k]])
                    else:
                        p_collection.set_verts([vtx_rot[k:k + num_vtx] for k in range(0, len(vtx_rot), num_vtx)])

                    h_lines[ORIENTATION][1].set_xdata([x_min, x_max])
                    h_lines[ORIENTATION][1].set_ydata([0, 0])
                    h_lines[ORIENTATION][1].set_3d_properties([0, 0])
                    h_lines[ORIENTATION][2].set_xdata([0, 0])
                    h_lines[ORIENTATION][2].set_ydata([y_min, y_max])
                    h_lines[ORIENTATION][2].set_3d_properties([0, 0])
                    h_lines[ORIENTATION][3].set_xdata([0, 0])
                    h_lines[ORIENTATION][3].set_ydata([0, 0])
                    h_lines[ORIENTATION][3].set_3d_properties([z_min, z_max])

                    need_redraw = True

                elif i == R:
                    h_lines[R][0].set_xdata(data[R][::3])
                    h_lines[R][0].set_ydata(data[R][1::3])
                    h_lines[R][0].set_3d_properties(data[R][2::3])

                    x_min = min(data[R][::3])
                    x_max = max(data[R][::3])
                    x_center = 0.5 * (x_max + x_min)

                    y_min = min(data[R][1::3])
                    y_max = max(data[R][1::3])
                    y_center = 0.5 * (y_max + y_min)

                    z_min = min(data[R][2::3])
                    z_max = max(data[R][2::3])
                    z_center = 0.5 * (z_max + z_min)

                    total_min = min([x_min, y_min, z_min])
                    total_max = max([x_max, y_max, z_max])
                    half_span = 0.5 * (total_max - total_min)

                    h_plots[R].set_xlim3d(x_center - half_span, x_center + half_span)
                    h_plots[R].set_ylim3d(y_center - half_span, y_center + half_span)
                    h_plots[R].set_zlim3d(z_center - half_span, z_center + half_span)
                    need_redraw = True
                else:
                    for k in range(elems[i]):
                        h_lines[i][k].set_data(data[0], data[i][k::elems[i]])

                    h_fig.canvas.restore_region(canvas[i])
                    for artist in h_lines[i]:
                        h_plots[i].draw_artist(artist)
                    h_fig.canvas.blit(h_plots[i].bbox)

                    old_lim = h_plots[i].get_xlim(), h_plots[i].get_ylim()
                    h_plots[i].relim()
                    h_plots[i].autoscale_view()
                    if (h_plots[i].get_xlim(), h_plots[i].get_ylim()) != old_lim:
                        need_redraw = True

            if need_redraw:
                h_fig.canvas.draw()

            h_fig.canvas.flush_events()

# Any of the drawing commands may fail because the window
# has closed, so swallow those exceptions. We're closing anyway.
except (AttributeError, TclError) as e:
    onclose(None)
