#!/usr/bin/env python3.4

####################################################################
#   KPS_Vis.py
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
#   Realtime visualization of any combination of the eight parameters
#   selectable below. Just comment out the ones you don't want.
#
#   Visualization begins in realtime mode, constantly checking
#   the outfiles for updates. When no changes occur in the files
#   for some time, KPS_Vis switches to FINAL mode, adjusting the axes
#   and allowing user interaction.
#

PLOT_R, PLOT_V, PLOT_Q, PLOT_Q_ORB, PLOT_W, PLOT_V_B, PLOT_E, PLOT_ALT, PLOT_B_STAR, PLOT_ORIENTATION = [False] * 10
#####################################################################
##################  USER CONFIGURABLE PARAMETERS  ###################
#####################################################################

font_size = 14
line_width = 1.5
maximize_plot = True

# for live orientation plotter
wireframe = False
poly_color = 'blue'
num_vtx = 4
poly_file = 'poly.kps'
face_alpha = 0.8
axis_alpha = 0.5

# Don't change these to False; simply comment
# out the ones you don't want.

#PLOT_R = True
#PLOT_V = True
#PLOT_Q = True
#PLOT_Q_ORB = True
#PLOT_ALT = True
#PLOT_B_STAR = True
#PLOT_W = True
#PLOT_V_B = True
#PLOT_E = True
PLOT_ORIENTATION = True

#####################################################################
#####################################################################
#####################################################################





import sys
from math import ceil
import matplotlib.pyplot as plt
import numpy as np
from tkinter import TclError

if PLOT_R or PLOT_ORIENTATION:
    from mpl_toolkits.mplot3d import Axes3D

if PLOT_ORIENTATION:
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

def initPlots():
    global vtx, r_plt, v_plt, v_b_plt, q_plt, q_orb_plt, w_plt, e_plt, alt_plt, b_star_plt, or_plt, or_0, or_1, or_2, or_3, r_0, r_1, q_0, q_1, q_2, q_3, q_orb_0, q_orb_1, q_orb_2, q_orb_3, w_x, w_y, w_z, v_x, v_y, v_z, v_b_x, v_b_y, v_b_z, e_0, alt_0, b_star_0, w_canvas, v_canvas, v_b_canvas, e_canvas, q_canvas, q_orb_canvas, alt_canvas, b_star_canvas
    cur_plot = 1

    if PLOT_R:
        r_plt = fig.add_subplot(num_plots, 1, cur_plot, aspect='equal', projection='3d')
        cur_plot += 1
        r_plt.set_title('Position', fontweight="bold")
        r_plt.set_xlabel(r'$x_{ECI}  [m]$')
        r_plt.set_ylabel(r'$y_{ECI}  [m]$')
        r_plt.set_zlabel(r'$z_{ECI}  [m]$')
        r_0 = r_plt.plot([], [], [], color='g', linewidth=line_width)[0]
        r_1 = r_plt.plot([0], [0], [0], 'bo', markersize=40, alpha=0.90)[0]

    if PLOT_ORIENTATION:
        or_plt = fig.add_subplot(num_plots, 1, cur_plot, aspect='equal', projection='3d')
        cur_plot += 1
        or_plt.set_title('Live Orientation', fontweight="bold")
        or_plt.set_xlabel(r'$x_{orbital}  [m]$')
        or_plt.set_ylabel(r'$y_{orbital}  [m]$')
        or_plt.set_zlabel(r'$z_{orbital}  [m]$')
        if wireframe:
            for i in range(num_poly):
                or_0.append(or_plt.plot([], [], [], color=poly_color, linewidth=line_width)[0])
        else:
            or_0 = Poly3DCollection(or_0, alpha=face_alpha, edgecolor='k', color=poly_color)
            or_plt.add_collection3d(or_0)

        or_1 = or_plt.plot([], [], [], color='k', alpha=axis_alpha)[0]
        or_2 = or_plt.plot([], [], [], color='k', alpha=axis_alpha)[0]
        or_3 = or_plt.plot([], [], [], color='k', alpha=axis_alpha)[0]

    if PLOT_V:
        v_plt = fig.add_subplot(num_plots, 1, cur_plot)
        cur_plot += 1
        v_plt.set_title('Velocity', fontweight="bold")
        v_plt.set_ylabel('[m/s]')
        v_x = v_plt.plot([], [], color='r', label='$v_x$', linewidth=line_width)[0]
        v_y = v_plt.plot([], [], color='g', label='$v_y$', linewidth=line_width)[0]
        v_z = v_plt.plot([], [], color='b', label='$v_z$', linewidth=line_width)[0]
        v_plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    if PLOT_Q:
        q_plt = fig.add_subplot(num_plots, 1, cur_plot)
        cur_plot += 1
        q_plt.set_title('Quaternion to ECI Frame', fontweight="bold")
        q_plt.set_ylabel('Component')
        q_0 = q_plt.plot([], [], color='r', label='$q_0$', linewidth=line_width)[0]
        q_1 = q_plt.plot([], [], color='g', label='$q_1$', linewidth=line_width)[0]
        q_2 = q_plt.plot([], [], color='b', label='$q_2$', linewidth=line_width)[0]
        q_3 = q_plt.plot([], [], color='k', label='$q_3$', linewidth=line_width)[0]
        q_plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    if PLOT_Q_ORB:
        q_orb_plt = fig.add_subplot(num_plots, 1, cur_plot)
        cur_plot += 1
        q_orb_plt.set_title('Quaternion to Orbital Frame', fontweight="bold")
        q_orb_plt.set_ylabel('Component')
        q_orb_0 = q_orb_plt.plot([], [], color='r', label='$q_0$', linewidth=line_width)[0]
        q_orb_1 = q_orb_plt.plot([], [], color='g', label='$q_1$', linewidth=line_width)[0]
        q_orb_2 = q_orb_plt.plot([], [], color='b', label='$q_2$', linewidth=line_width)[0]
        q_orb_3 = q_orb_plt.plot([], [], color='k', label='$q_3$', linewidth=line_width)[0]
        q_orb_plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    if PLOT_W:
        w_plt = fig.add_subplot(num_plots, 1, cur_plot)
        cur_plot += 1
        w_plt.set_title('Angular Velocity', fontweight="bold")
        w_plt.set_ylabel('$[s^{-1}]$')
        w_x = w_plt.plot([], [], color='r', label='$\omega_x$', linewidth=line_width)[0]
        w_y = w_plt.plot([], [], color='g', label='$\omega_y$', linewidth=line_width)[0]
        w_z = w_plt.plot([], [], color='b', label='$\omega_z$', linewidth=line_width)[0]
        w_plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    if PLOT_V_B:
        v_b_plt = fig.add_subplot(num_plots, 1, cur_plot)
        cur_plot += 1
        v_b_plt.set_title('Velocity in Body Frame', fontweight="bold")
        v_b_plt.set_ylabel('[m/s]')
        v_b_x = v_b_plt.plot([], [], color='r', label='$v_x$', linewidth=line_width)[0]
        v_b_y = v_b_plt.plot([], [], color='g', label='$v_y$', linewidth=line_width)[0]
        v_b_z = v_b_plt.plot([], [], color='b', label='$v_z$', linewidth=line_width)[0]
        v_b_plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    if PLOT_E:
        e_plt = fig.add_subplot(num_plots, 1, cur_plot)
        cur_plot += 1
        e_plt.set_title('Pointing Error', fontweight="bold")
        e_plt.set_ylabel('[' + chr(176) + ']')
        e_0 = e_plt.plot([], [], color='r', label='$\epsilon_p$', linewidth=line_width)[0]
        e_plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        e_plt.set_ylim(0, 180)

    if PLOT_ALT:
        alt_plt = fig.add_subplot(num_plots, 1, cur_plot)
        cur_plot += 1
        alt_plt.set_title('Altitude', fontweight="bold")
        alt_plt.set_ylabel('[m]')
        alt_0 = alt_plt.plot([], [], color='g', label='alt', linewidth=line_width)[0]
        alt_plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    if PLOT_B_STAR:
        b_star_plt = fig.add_subplot(num_plots, 1, cur_plot)
        cur_plot += 1
        b_star_plt.set_title('B$^*$ Term', fontweight="bold")
        b_star_plt.set_ylabel('[$m^{-1}$]')
        b_star_0 = b_star_plt.plot([], [], color='b', label='$B^*$', linewidth=line_width)[0]
        b_star_plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    if (PLOT_R and num_plots == 1) or (PLOT_ORIENTATION and num_plots == 1) or (PLOT_R and PLOT_ORIENTATION and num_plots == 2):
        fig.tight_layout()
    else:
        plt.xlabel('Time [s]', fontsize=1.4*font_size)

    fig.canvas.draw()

    # store bbox in order to blit and restore later (for speed)
    if PLOT_Q:
        q_canvas = fig.canvas.copy_from_bbox(q_plt.bbox)

    if PLOT_Q_ORB:
        q_orb_canvas = fig.canvas.copy_from_bbox(q_orb_plt.bbox)

    if PLOT_W:
        w_canvas = fig.canvas.copy_from_bbox(w_plt.bbox)

    if PLOT_V:
        v_canvas = fig.canvas.copy_from_bbox(v_plt.bbox)

    if PLOT_V_B:
        v_b_canvas = fig.canvas.copy_from_bbox(v_b_plt.bbox)

    if PLOT_E:
        e_canvas = fig.canvas.copy_from_bbox(e_plt.bbox)

    if PLOT_ALT:
        alt_canvas = fig.canvas.copy_from_bbox(alt_plt.bbox)

    if PLOT_B_STAR:
        b_star_canvas = fig.canvas.copy_from_bbox(b_star_plt.bbox)

num_plots = PLOT_R + PLOT_V + PLOT_Q + PLOT_Q_ORB + PLOT_ORIENTATION + PLOT_W + PLOT_V_B + PLOT_E + PLOT_ALT + PLOT_B_STAR

# read all data available so far

try:
    f_t = open('t.bin', 'rb')
except:
    print('ERROR: failed to open file. Aborting.')
    sys.exit(-1)
# discard last element because it may be incomplete, as it's being written in realtime
t = np.fromfile(f_t)[:-1]

if PLOT_R:
    try:
        f_r = open('r.bin', 'rb')
    except:
        print('ERROR: failed to open file. Aborting.')
        sys.exit(-1)
    r = np.fromfile(f_r)
    r = r[:3*ceil(r.size/3) - 2]

if PLOT_V:
    try:
        f_v = open('v.bin', 'rb')
    except:
        print('ERROR: failed to open file. Aborting.')
        sys.exit(-1)
    v = np.fromfile(f_v)
    v = v[:3*ceil(v.size/3) - 2]

if PLOT_V_B:
    try:
        f_v_b = open('v_body.bin', 'rb')
    except:
        print('ERROR: failed to open file. Aborting.')
        sys.exit(-1)
    v_b = np.fromfile(f_v_b)
    v_b = v_b[:3*ceil(len(v_b)/3) - 2]

if PLOT_E:
    try:
        f_e = open('p_e.bin', 'rb')
    except:
        print('ERROR: failed to open file. Aborting.')
        sys.exit(-1)
    e = np.fromfile(f_e)[:-1]

if PLOT_ALT:
    try:
        f_alt = open('alt.bin', 'rb')
    except:
        print('ERROR: failed to open file. Aborting.')
        sys.exit(-1)
    alt = np.fromfile(f_alt)[:-1]

if PLOT_B_STAR:
    try:
        f_b_star = open('b_star.bin', 'rb')
    except:
        print('ERROR: failed to open file. Aborting.')
        sys.exit(-1)
    b_star = np.fromfile(f_b_star)[:-1]

if PLOT_Q:
    try:
        f_q = open('q.bin', 'rb')
    except:
        print('ERROR: failed to open file. Aborting.')
        sys.exit(-1)
    q = np.fromfile(f_q)
    q = q[:4*ceil(q.size/4) - 3]

if PLOT_Q_ORB or PLOT_ORIENTATION:
    try:
        f_q_orb = open('q_orb.bin', 'rb')
    except:
        print('ERROR: failed to open file. Aborting.')
        sys.exit(-1)
    q_orb = np.fromfile(f_q_orb)
    q_orb = q_orb[:4*ceil(q_orb.size/4) - 3]

if PLOT_W:
    try:
        f_w = open('w.bin', 'rb')
    except:
        print('ERROR: failed to open file. Aborting.')
        sys.exit(-1)
    w = np.fromfile(f_w)
    w = w[:3*ceil(w.size/3) - 2]

if PLOT_ORIENTATION:
    try:
        f_orb = open(poly_file)
    except:
        print('ERROR: failed to open polygon file. Aborting.')
        sys.exit(-1)
    vtx = np.array([[float(x) for x in line.rstrip('\n').split(', ')] for line in f_orb if len(line) > 4 and line[0] != '#'])
    num_poly = len(vtx) // num_vtx;
    or_0 = []
    f_orb.close()

plt.rcParams.update({'font.size': font_size})

fig = plt.figure('KPS - REALTIME')

r_plt, v_plt, v_b_plt, q_plt, q_orb_plt, w_plt, e_plt, alt_plt, b_star_plt, or_plt = [None] * 10

if maximize_plot:
    maximizePlot()

plt.show(False)

initPlots()

ticks_without_change = 0

close_please = False

def onresize(event):
    fig.clf()
    initPlots()

def onclose(event):
    global close_please
    close_please = True

# force redraw on resize
cid = fig.canvas.mpl_connect('resize_event', onresize)

# graceful close (matplotlib makes this difficult)
fig.canvas.mpl_connect('close_event', onclose)

while True:
    if close_please:
        break
    try:
        f_t.seek(8*t.size)
        t_in = np.fromfile(f_t)[:-1]
        if t_in.size:
            ticks_without_change = 0
        else:
            ticks_without_change += 1
            if ticks_without_change > 20:
                # no change for a while
                # switch to Final mode
                fig.canvas.set_window_title('KPS - Final')
                fig.canvas.mpl_disconnect(cid)
                if PLOT_V:
                    v_plt.relim()
                    v_plt.autoscale_view(True, True, False)
                if PLOT_V_B:
                    v_b_plt.relim()
                    v_b_plt.autoscale_view(True, True, False)
                if PLOT_W:
                    w_plt.relim()
                    w_plt.autoscale_view(True, True, False)
                if PLOT_Q:
                    q_plt.relim()
                    q_plt.autoscale_view(True, True, False)
                if PLOT_Q_ORB:
                    q_orb_plt.relim()
                    q_orb_plt.autoscale_view(True, True, False)
                if PLOT_E:
                    e_plt.relim()
                    e_plt.autoscale_view(True, True, False)
                if PLOT_ALT:
                    alt_plt.relim()
                    alt_plt.autoscale_view(True, True, False)
                if PLOT_B_STAR:
                    b_star_plt.relim()
                    b_star_plt.autoscale_view(True, True, False)
                plt.show()
        t = np.append(t, t_in)

        if PLOT_R:
            f_r.seek(8*r.size)
            r_in = np.fromfile(f_r)
            r = np.append(r, r_in[:3*ceil(r_in.size/3) - 2])

        if PLOT_V:
            fig.canvas.restore_region(v_canvas)
            v_plt.draw_artist(v_x)
            v_plt.draw_artist(v_y)
            v_plt.draw_artist(v_z)
            fig.canvas.blit(v_plt.bbox)
            f_v.seek(8*v.size)
            v_in = np.fromfile(f_v)
            v = np.append(v, v_in[:3*ceil(v_in.size/3) - 2])

        if PLOT_V_B:
            fig.canvas.restore_region(v_b_canvas)
            v_b_plt.draw_artist(v_b_x)
            v_b_plt.draw_artist(v_b_y)
            v_b_plt.draw_artist(v_b_z)
            fig.canvas.blit(v_b_plt.bbox)
            f_v_b.seek(8*len(v_b))
            v_b_in = np.fromfile(f_v_b)
            v_b = np.append(v_b, v_b_in[:3*ceil(v_b_in.size/3) - 2])

        if PLOT_E:
            fig.canvas.restore_region(e_canvas)
            e_plt.draw_artist(e_0)
            fig.canvas.blit(e_plt.bbox)
            f_e.seek(8*e.size)
            e = np.append(e, np.fromfile(f_e)[:-1])

        if PLOT_ALT:
            fig.canvas.restore_region(alt_canvas)
            alt_plt.draw_artist(alt_0)
            fig.canvas.blit(alt_plt.bbox)
            f_alt.seek(8*alt.size)
            alt = np.append(alt, np.fromfile(f_alt)[:-1])

        if PLOT_B_STAR:
            fig.canvas.restore_region(b_star_canvas)
            b_star_plt.draw_artist(b_star_0)
            fig.canvas.blit(b_star_plt.bbox)
            f_b_star.seek(8*b_star.size)
            b_star = np.append(b_star, np.fromfile(f_b_star)[:-1])

        if PLOT_W:
            fig.canvas.restore_region(w_canvas)
            w_plt.draw_artist(w_x)
            w_plt.draw_artist(w_y)
            w_plt.draw_artist(w_z)
            fig.canvas.blit(w_plt.bbox)
            f_w.seek(8*w.size)
            w_in = np.fromfile(f_w)
            w = np.append(w, w_in[:3*ceil(w_in.size/3) - 2])

        if PLOT_Q:
            fig.canvas.restore_region(q_canvas)
            q_plt.draw_artist(q_0)
            q_plt.draw_artist(q_1)
            q_plt.draw_artist(q_2)
            q_plt.draw_artist(q_3)
            fig.canvas.blit(q_plt.bbox)
            f_q.seek(8*q.size)
            q_in = np.fromfile(f_q)
            q = np.append(q, q_in[:4*ceil(q_in.size/4) - 3])

        if PLOT_Q_ORB:
            fig.canvas.restore_region(q_orb_canvas)
            q_orb_plt.draw_artist(q_orb_0)
            q_orb_plt.draw_artist(q_orb_1)
            q_orb_plt.draw_artist(q_orb_2)
            q_orb_plt.draw_artist(q_orb_3)
            fig.canvas.blit(q_orb_plt.bbox)

        if PLOT_Q_ORB or PLOT_ORIENTATION:
            f_q_orb.seek(8*q_orb.size)
            q_orb_in = np.fromfile(f_q_orb)
            q_orb = np.append(q_orb, q_orb_in[:4*ceil(q_orb_in.size/4) - 3])

        # determine length of shortest file
        # (all files will be trimmed to this length for plotting)
        length = t.size
        if PLOT_R:
            if r.size//3 < length:
                length = r.size//3;
        if PLOT_V:
            if v.size//3 < length:
                length = v.size//3;
        if PLOT_Q:
            if q.size//4 < length:
                length = q.size//4;
        if PLOT_Q_ORB or PLOT_ORIENTATION:
            if q_orb.size//4 < length:
                length = q_orb.size//4;
        if PLOT_W:
            if w.size//3 < length:
                length = w.size//3;
        if PLOT_V_B:
            if len(v_b)//3 < length:
                length = len(v_b)//3;
        if PLOT_E:
            if e.size < length:
                length = e.size;

        t = t[:length]        

        need_redraw = False

        if PLOT_R:
            r = r[:3*length]
            r_0.set_xdata(r[::3])
            r_0.set_ydata(r[1::3])
            r_0.set_3d_properties(r[2::3])

            x_min = min(r[::3])
            x_max = max(r[::3])
            x_center = 0.5 * (x_max + x_min)

            y_min = min(r[1::3])
            y_max = max(r[1::3])
            y_center = 0.5 * (y_max + y_min)

            z_min = min(r[2::3])
            z_max = max(r[2::3])
            z_center = 0.5 * (z_max + z_min)

            total_min = min([x_min, y_min, z_min])
            total_max = max([x_max, y_max, z_max])
            half_span = 0.5 * (total_max - total_min)

            r_plt.set_xlim3d(x_center - half_span, x_center + half_span)
            r_plt.set_ylim3d(y_center - half_span, y_center + half_span)
            r_plt.set_zlim3d(z_center - half_span, z_center + half_span)

            need_redraw = True

        if PLOT_ORIENTATION:      
            q_vec = q_orb[-3:]
            vtx_rot = [v + cross(2*q_vec,cross(q_vec, v)+q_orb[-4]*v) for v in vtx]

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

            or_plt.set_xlim3d(x_center - half_span, x_center + half_span)
            or_plt.set_ylim3d(y_center - half_span, y_center + half_span)
            or_plt.set_zlim3d(z_center - half_span, z_center + half_span)

            if wireframe:
                for i in range(0, len(z), num_vtx):
                    or_0[i//num_vtx].set_xdata(x[i:i + num_vtx] + [x[i]])
                    or_0[i//num_vtx].set_ydata(y[i:i + num_vtx] + [y[i]])
                    or_0[i//num_vtx].set_3d_properties(z[i:i + num_vtx] + [z[i]])
            else:
                or_0.set_verts([vtx_rot[i:i + num_vtx] for i in range(0, len(vtx_rot), num_vtx)])

            or_1.set_xdata([x_min, x_max])
            or_1.set_ydata([0, 0])
            or_1.set_3d_properties([0, 0])
            or_2.set_xdata([0, 0])
            or_2.set_ydata([y_min, y_max])
            or_2.set_3d_properties([0, 0])
            or_3.set_xdata([0, 0])
            or_3.set_ydata([0, 0])
            or_3.set_3d_properties([z_min, z_max])

            need_redraw = True

        if PLOT_V:
            v = v[:3*length]
            v_x.set_data(t, v[::3])
            v_y.set_data(t, v[1::3])
            v_z.set_data(t, v[2::3])
            old_lim = v_plt.get_xlim(), v_plt.get_ylim()
            v_plt.relim()
            v_plt.autoscale_view()
            if (v_plt.get_xlim(), v_plt.get_ylim()) != old_lim:
                need_redraw = True

        if PLOT_Q:
            q = q[:4*length]
            q_0.set_data(t, q[::4])
            q_1.set_data(t, q[1::4])
            q_2.set_data(t, q[2::4])
            q_3.set_data(t, q[3::4])
            old_lim = q_plt.get_xlim(), q_plt.get_ylim()
            q_plt.relim()
            q_plt.autoscale_view()
            if (q_plt.get_xlim(), q_plt.get_ylim()) != old_lim:
                need_redraw = True

        if PLOT_Q_ORB:
            q_orb = q_orb[:4*length]
            q_orb_0.set_data(t, q_orb[::4])
            q_orb_1.set_data(t, q_orb[1::4])
            q_orb_2.set_data(t, q_orb[2::4])
            q_orb_3.set_data(t, q_orb[3::4])
            old_lim = q_orb_plt.get_xlim(), q_orb_plt.get_ylim()
            q_orb_plt.relim()
            q_orb_plt.autoscale_view()
            if (q_orb_plt.get_xlim(), q_orb_plt.get_ylim()) != old_lim:
                need_redraw = True

        if PLOT_W:
            w = w[:3*length]
            w_x.set_data(t, w[::3])
            w_y.set_data(t, w[1::3])
            w_z.set_data(t, w[2::3])
            old_lim = w_plt.get_xlim(), w_plt.get_ylim()
            w_plt.relim()
            w_plt.autoscale_view()
            if (w_plt.get_xlim(), w_plt.get_ylim()) != old_lim:
                need_redraw = True

        if PLOT_V_B:
            v_b = v_b[:3*length]
            v_b_x.set_data(t, v_b[::3])
            v_b_y.set_data(t, v_b[1::3])
            v_b_z.set_data(t, v_b[2::3])
            old_lim = v_b_plt.get_xlim(), v_b_plt.get_ylim()
            v_b_plt.relim()
            v_b_plt.autoscale_view()
            if (v_b_plt.get_xlim(), v_b_plt.get_ylim()) != old_lim:
                need_redraw = True

        if PLOT_E:
            e = e[:length]
            e_0.set_data(t, e)
            old_lim = e_plt.get_xlim(), e_plt.get_ylim()
            e_plt.relim()
            e_plt.autoscale_view()
            if (e_plt.get_xlim(), e_plt.get_ylim()) != old_lim:
                need_redraw = True

        if PLOT_ALT:
            alt = alt[:length]
            alt_0.set_data(t, alt)
            old_lim = alt_plt.get_xlim(), alt_plt.get_ylim()
            alt_plt.relim()
            alt_plt.autoscale_view()
            if (alt_plt.get_xlim(), alt_plt.get_ylim()) != old_lim:
                need_redraw = True

        if PLOT_B_STAR:
            b_star = b_star[:length]
            b_star_0.set_data(t, b_star)
            old_lim = b_star_plt.get_xlim(), b_star_plt.get_ylim()
            b_star_plt.relim()
            b_star_plt.autoscale_view()
            if (b_star_plt.get_xlim(), b_star_plt.get_ylim()) != old_lim:
                need_redraw = True

        if need_redraw:
            fig.canvas.draw()

        fig.canvas.flush_events()

        # Tcl is pretty buggy. Realtime data makes it crash on the way down,
        # so we swallow the exception. We're closing anyway.
    except (AttributeError, TclError) as e:
        onclose(None)
