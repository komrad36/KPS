%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   KPS_Vis.m
%   KPS
%	
%	Author: Kareem Omar
%	kareem.omar@uah.edu
%	https://github.com/komrad36
%
%	Last updated Feb 12, 2016
%   This application is entirely my own work.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Realtime visualization of any combination of the eight parameters
%   selectable below. Just comment out the ones you don't want.
%
%   Visualization begins in realtime mode, constantly checking
%   the outfiles for updates. When no changes occur in the files
%   for some time, KPS_Vis switches to FINAL mode, adjusting the axes
%   and allowing user interaction.
%

function KPS_Vis
 global NUM_PLOTS AXIS_LINE_WIDTH FONT_SIZE LINE_WIDTH t r v v_b p_e q q_orb w alt b_star PLOT_R PLOT_V PLOT_Q PLOT_Q_ORB PLOT_W PLOT_V_B PLOT_E PLOT_ALT PLOT_B_STAR PLOT_ORIENTATION v_plt q_plt q_orb_plt w_plt v_b_plt e_plt alt_plt b_star_plt P num_vtx wireframe color_if_not_gradient gradient_if_not_wireframe face_alpha at_origin_axis_line_widths
initBools
%% User configurables

FONT_SIZE = 16;
LEGEND_FONT_SIZE = 14;
LINE_WIDTH = 1.2;
AXIS_LINE_WIDTH = 2;
MAXIMIZE_PLOT = true;

% for live orientation plotter
wireframe = false;
at_origin_axis_line_widths = 1;
gradient_if_not_wireframe = true;
color_if_not_gradient = 'blue';
num_vtx = 4;
poly_file = 'poly.kps';

% for MATLAB only, won't work in Octave
face_alpha = 0.8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('OCTAVE_VERSION', 'builtin')
  face_alpha = 1.0;
  LEGEND_FONT_SIZE = 8;
  graphics_toolkit('fltk')
end %if

%%% Don't change these to false; simply comment
%%% out the ones you don't want

% PLOT_R = true;
% PLOT_V = true;
% PLOT_Q = true;
% PLOT_Q_ORB = true;
% PLOT_ALT = true;
% PLOT_B_STAR = true;
% PLOT_W = true;
% PLOT_V_B = true;
% PLOT_E = true;
PLOT_ORIENTATION = true;


%% Execution

NUM_PLOTS = PLOT_R + PLOT_V + PLOT_Q + PLOT_Q_ORB + PLOT_W + PLOT_V_B + PLOT_E + PLOT_ALT + PLOT_B_STAR + PLOT_ORIENTATION;

f_t = fopen('t.bin', 'r');
if (f_t == -1)
    error('Failed to open file.')
end %if
t = fread(f_t, Inf, 'double');
t = t(1:end-1);
% discard last element because it may be incomplete, as it's being written in realtime

if PLOT_R
    f_r = fopen('r.bin', 'r');
    if (f_r == -1)
        error('Failed to open file.')
    end %if
    r = fread(f_r, Inf, 'double');
    r = r(1:3*ceil(numel(r)/3) - 3);
end %if

if PLOT_V
    f_v = fopen('v.bin', 'r');
    if (f_v == -1)
        error('Failed to open file.')
    end %if
    v = fread(f_v, Inf, 'double');
    v = v(1:3*ceil(numel(v)/3) - 3);
end %if

if PLOT_V_B
    f_v_b = fopen('v_body.bin', 'r');
    if (f_v_b == -1)
        error('Failed to open file.')
    end %if
    v_b = fread(f_v_b, Inf, 'double');
    v_b = v_b(1:3*ceil(numel(v_b)/3) - 3);
end %if

if PLOT_E
    f_p_e = fopen('p_e.bin', 'r');
    if (f_p_e == -1)
        error('Failed to open file.')
    end %if
    p_e = fread(f_p_e, Inf, 'double');
    p_e = p_e(1:end-1);
end %if

if PLOT_ALT
    f_alt = fopen('alt.bin', 'r');
    if (f_alt == -1)
        error('Failed to open file.')
    end %if
    alt = fread(f_alt, Inf, 'double');
    alt = alt(1:end-1);
end %if

if PLOT_B_STAR
    f_b_star = fopen('b_star.bin', 'r');
    if (f_b_star == -1)
        error('Failed to open file.')
    end %if
    b_star = fread(f_b_star, Inf, 'double');
    b_star = b_star(1:end-1);
end %if

if PLOT_Q
    f_q = fopen('q.bin', 'r');
    if (f_q == -1)
        error('Failed to open file.')
    end %if
    q = fread(f_q, Inf, 'double');
    q = q(1:4*ceil(numel(q)/4) - 4);
end %if

if PLOT_Q_ORB || PLOT_ORIENTATION
    f_q_orb = fopen('q_orb.bin', 'r');
    if (f_q_orb == -1)
        error('Failed to open file.')
    end %if
    q_orb = fread(f_q_orb, Inf, 'double');
    q_orb = q_orb(1:4*ceil(numel(q_orb)/4) - 4);
end %if

if PLOT_W
    f_w = fopen('w.bin', 'r');
    if (f_w == -1)
        error('Failed to open file.')
    end %if
    w = fread(f_w, Inf, 'double');
    w = w(1:3*ceil(numel(w)/3) - 3);
end %if

if PLOT_ORIENTATION
    fid = fopen(poly_file, 'r');
    if fid == -1
       error('Unable to open polygon file. Aborting.')  
    end %if
    
    % load polygons
    % this could be textscan() in MATLAB but Octave
    % doesn't yet support '%[]' (skip characters) arguments
    % so it's handrolled for now
    P = zeros(3, 0);
    i = 1;
    while true
      ln = fgetl(fid);
      if ln == -1, break, end

      if numel(ln) == 0 || ln(1) == '#', continue, end
      elems = str2double(strsplit(ln, ','));
      if numel(elems) ~= 3
        error('Error parsing file!')
      end
      P(:, i) = elems;
      i = i + 1;
    end %while

    fclose(fid);
    clear fid
end %if

h = figure('Name', 'KPS - REALTIME', 'NumberTitle', 'off'); clf
if MAXIMIZE_PLOT
    set(h, 'units','normalized','outerposition',[0 0 1 1])
end %if

initPlots
updatePlots

if PLOT_Q, l=legend(q_plt, 'q_0', 'q_1', 'q_2', 'q_3', 'location', 'eastoutside'); set(l, 'FontSize', LEGEND_FONT_SIZE), end
if PLOT_Q_ORB, l=legend(q_orb_plt, 'q_0', 'q_1', 'q_2', 'q_3', 'location', 'eastoutside'); set(l, 'FontSize', LEGEND_FONT_SIZE), end
if PLOT_W, l=legend(w_plt, '\omega_x', '\omega_y', '\omega_z', 'location', 'eastoutside'); set(l, 'FontSize', LEGEND_FONT_SIZE), end
if PLOT_V, l=legend(v_plt, 'v_x', 'v_y', 'v_z', 'location', 'eastoutside'); set(l, 'FontSize', LEGEND_FONT_SIZE), end
if PLOT_V_B, l=legend(v_b_plt, 'v_x', 'v_y', 'v_z', 'location', 'eastoutside'); set(l, 'FontSize', LEGEND_FONT_SIZE), end
if PLOT_E, l=legend(e_plt, '\epsilon_p', 'location', 'eastoutside'); set(l, 'FontSize', LEGEND_FONT_SIZE), end
if PLOT_ALT, l=legend(alt_plt, 'alt', 'location', 'eastoutside'); set(l, 'FontSize', LEGEND_FONT_SIZE), end
if PLOT_B_STAR, l=legend(b_star_plt, 'B*', 'location', 'eastoutside'); set(l, 'FontSize', LEGEND_FONT_SIZE), end

ticks_without_change = 0;

while true
    if ~ishghandle(h), break, end
    drawnow
    
    fseek(f_t, 8*numel(t), 'bof');
    t_in = fread(f_t, Inf, 'double');
    t_in(end) = [];
    if numel(t_in)
        ticks_without_change = 0;
    else
        ticks_without_change = ticks_without_change + 1;
        if ticks_without_change > 20
            % no change for a while
            % switch to Final mode
            fixAxes
            set(h, 'Name', 'KPS - Final')
            break
        end %if
    end %if
    t = [t; t_in];

    if PLOT_R
        fseek(f_r, 8*numel(r), 'bof');
        r_in = fread(f_r, Inf, 'double');
        r_in(3*ceil(numel(r_in)/3) - 2:end) = [];
        r = [r; r_in];
    end %if
    
    if PLOT_V
        fseek(f_v, 8*numel(v), 'bof');
        v_in = fread(f_v, Inf, 'double');
        v_in(3*ceil(numel(v_in)/3) - 2:end) = [];
        v = [v; v_in];
    end %if
    
    if PLOT_V_B
        fseek(f_v_b, 8*numel(v_b), 'bof');
        v_b_in = fread(f_v_b, Inf, 'double');
        v_b_in(3*ceil(numel(v_b_in)/3) - 2:end) = [];
        v_b = [v_b; v_b_in];
    end %if
    
    if PLOT_E
        fseek(f_p_e, 8*numel(p_e), 'bof');
        p_e_in = fread(f_p_e, Inf, 'double');
        p_e_in(end) = [];
        p_e = [p_e; p_e_in];
    end %if
    
    if PLOT_ALT
        fseek(f_alt, 8*numel(alt), 'bof');
        alt_in = fread(f_alt, Inf, 'double');
        alt_in(end) = [];
        alt = [alt; alt_in];
    end %if
    
    if PLOT_B_STAR
        fseek(f_b_star, 8*numel(b_star), 'bof');
        b_star_in = fread(f_b_star, Inf, 'double');
        b_star_in(end) = [];
        b_star = [b_star; b_star_in];
    end %if
    
    if PLOT_Q
        fseek(f_q, 8*numel(q), 'bof');
        q_in = fread(f_q, Inf, 'double');
        q_in(4*ceil(numel(q_in)/4) - 3:end) = [];
        q = [q; q_in];
    end %if
    
    if PLOT_Q_ORB || PLOT_ORIENTATION
        fseek(f_q_orb, 8*numel(q_orb), 'bof');
        q_orb_in = fread(f_q_orb, Inf, 'double');
        q_orb_in(4*ceil(numel(q_orb_in)/4) - 3:end) = [];
        q_orb = [q_orb; q_orb_in];
    end %if
    
    if PLOT_W
        fseek(f_w, 8*numel(w), 'bof');
        w_in = fread(f_w, Inf, 'double');
        w_in(3*ceil(numel(w_in)/3) - 2:end) = [];
        w = [w; w_in];
    end %if
    
    if ~ishghandle(h), break, end
    updatePlots
    
end %while
fclose('all');
end %function

function initPlots
global NUM_PLOTS FONT_SIZE AXIS_LINE_WIDTH LINE_WIDTH PLOT_R PLOT_V PLOT_Q PLOT_Q_ORB PLOT_W PLOT_V_B PLOT_E PLOT_ALT PLOT_B_STAR PLOT_ORIENTATION r_plt v_plt q_plt q_orb_plt w_plt v_b_plt e_plt alt_plt b_star_plt or_plt
cur_plot = 1;

if PLOT_R
    r_plt = subplot(NUM_PLOTS, 1, cur_plot); hold on
    cur_plot = cur_plot + 1;
    axis equal
    title('Position', 'FontSize', FONT_SIZE)
    xlabel('x_{ECI}')
    ylabel('y_{ECI}')
    zlabel('z_{ECI}')
    set(gca, 'FontSize', FONT_SIZE, 'LineWidth', AXIS_LINE_WIDTH);
end %if

if PLOT_ORIENTATION
    or_plt = subplot(NUM_PLOTS, 1, cur_plot); hold on
    cur_plot = cur_plot + 1;
    axis equal
    title('Live Orientation', 'FontSize', FONT_SIZE)
    xlabel('x_{orbital}')
    ylabel('y_{orbital}')
    zlabel('z_{orbital}')
    set(gca, 'FontSize', FONT_SIZE, 'LineWidth', AXIS_LINE_WIDTH);
end %if

if PLOT_V
    v_plt = subplot(NUM_PLOTS, 1, cur_plot); hold on
    cur_plot = cur_plot + 1;
    title('Velocity', 'FontSize', FONT_SIZE)
    ylabel('[m/s]', 'FontSize', FONT_SIZE)
    set(gca, 'FontSize', FONT_SIZE, 'LineWidth', LINE_WIDTH);
end %if

if PLOT_Q
    q_plt = subplot(NUM_PLOTS, 1, cur_plot); hold on
    cur_plot = cur_plot + 1;
    title('Quaternion to ECI Frame', 'FontSize', FONT_SIZE)
    ylabel('Component', 'FontSize', FONT_SIZE)
    set(gca, 'FontSize', FONT_SIZE, 'LineWidth', LINE_WIDTH);
end %if

if PLOT_Q_ORB
    q_orb_plt = subplot(NUM_PLOTS, 1, cur_plot); hold on
    cur_plot = cur_plot + 1;
    title('Quaternion to Orbital Frame', 'FontSize', FONT_SIZE)
    ylabel('Component', 'FontSize', FONT_SIZE)
    set(gca, 'FontSize', FONT_SIZE, 'LineWidth', LINE_WIDTH);
end %if
    
if PLOT_W
    w_plt = subplot(NUM_PLOTS, 1, cur_plot); hold on
    cur_plot = cur_plot + 1;
    title('Angular Velocity', 'FontSize', FONT_SIZE)
    ylabel('[s^-^1]', 'FontSize', FONT_SIZE)
    set(gca, 'FontSize', FONT_SIZE, 'LineWidth', LINE_WIDTH);
end %if

if PLOT_V_B
    v_b_plt = subplot(NUM_PLOTS, 1, cur_plot); hold on
    cur_plot = cur_plot + 1;
    title('Velocity in Body Frame', 'FontSize', FONT_SIZE)
    ylabel('[m/s]', 'FontSize', FONT_SIZE)
    set(gca, 'FontSize', FONT_SIZE, 'LineWidth', LINE_WIDTH);
end %if

if PLOT_E
    e_plt = subplot(NUM_PLOTS, 1, cur_plot); hold on
    title('Pointing Error', 'FontSize', FONT_SIZE)
    cur_plot = cur_plot + 1;
    ylabel('[\circ]', 'FontSize', FONT_SIZE)
    ylim([0 180])
    set(gca, 'FontSize', FONT_SIZE, 'LineWidth', LINE_WIDTH);
end %if

if PLOT_ALT
    alt_plt = subplot(NUM_PLOTS, 1, cur_plot); hold on
    cur_plot = cur_plot + 1;
    title('Altitude', 'FontSize', FONT_SIZE)
    ylabel('[m]', 'FontSize', FONT_SIZE)
    set(gca, 'FontSize', FONT_SIZE, 'LineWidth', LINE_WIDTH);
end %if

if PLOT_B_STAR
    b_star_plt = subplot(NUM_PLOTS, 1, cur_plot); hold on
    title('B* Term', 'FontSize', FONT_SIZE)
    ylabel('m^-^1', 'FontSize', FONT_SIZE)
    set(gca, 'FontSize', FONT_SIZE, 'LineWidth', LINE_WIDTH);
end %if
    
if ~((PLOT_R && NUM_PLOTS == 1) || (PLOT_ORIENTATION && NUM_PLOTS == 1) || (PLOT_R && PLOT_ORIENTATION && NUM_PLOTS == 2))
    xlabel('Time [s]', 'FontSize', FONT_SIZE)
end %if
end %function

function updatePlots
global LINE_WIDTH t r v v_b p_e q q_orb w alt b_star PLOT_R PLOT_V PLOT_Q PLOT_Q_ORB PLOT_W PLOT_V_B PLOT_E PLOT_ALT PLOT_B_STAR PLOT_ORIENTATION r_plt v_plt q_plt q_orb_plt w_plt v_b_plt e_plt alt_plt b_star_plt or_plt P num_vtx wireframe color_if_not_gradient gradient_if_not_wireframe face_alpha at_origin_axis_line_widths

% determine length of shortest file
% (all files will be trimmed to this length for plotting)
len = numel(t);
if PLOT_R
    if numel(r)/3 < len, len = numel(r)/3; end
end %if
if PLOT_V
    if numel(v)/3 < len, len = numel(v)/3; end
end %if
if PLOT_Q
    if numel(q)/4 < len, len = numel(q)/4; end
end %if
if PLOT_Q_ORB || PLOT_ORIENTATION
    if numel(q_orb)/4 < len, len = numel(q_orb)/4; end
end %if
if PLOT_W
    if numel(w)/3 < len, len = numel(w)/3; end
end %if
if PLOT_V_B
    if numel(v_b)/3 < len, len = numel(v_b)/3; end
end %if
if PLOT_E
    if numel(p_e) < len, len = numel(p_e); end
end %if
if PLOT_ALT
    if numel(alt) < len, len = numel(alt); end
end %if
if PLOT_B_STAR
    if numel(b_star) < len, len = numel(b_star); end
end %if
  
t(len+1:end) = [];

% MATLAB supports addpoints() for even more speed
% but Octave doesn't have addpoints(). Nor does
% MATLAB prior to R2015a, so this is best for
% compatability. It's still really fast due to
% careful design. Only the axes (points) are
% cleared on each frame, and 'hold' is not used,
% so only one set of points goes on each graph.
if PLOT_R
    r(3*len+1:end) = [];
    cla(r_plt)
    plot3(r_plt, 0,0,0, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 35)
    plot3(r_plt, r(1:3:end), r(2:3:end), r(3:3:end), 'g', 'LineWidth', 2);
end %if

if PLOT_Q_ORB || PLOT_ORIENTATION
    q_orb(4*len+1:end) = [];
end %if

if PLOT_ORIENTATION
    cla(or_plt)
    P_rot = zeros(size(P));
    for i = 1:size(P, 2)
        P_rot(:, i) = P(:, i) + cross(2*q_orb(end-2:end), cross(q_orb(end-2:end), P(:, i))+q_orb(end-3)*P(:, i));
    end %for
    num_poly = size(P, 2) / num_vtx;
    x = zeros(num_vtx, num_poly);
    y = zeros(num_vtx, num_poly);
    z = zeros(num_vtx, num_poly);
    for i = 1:num_poly
        x(:, i) = P_rot(1, num_vtx*(i-1)+1:num_vtx*i);
        y(:, i) = P_rot(2, num_vtx*(i-1)+1:num_vtx*i);
        z(:, i) = P_rot(3, num_vtx*(i-1)+1:num_vtx*i);
    end %if
    if wireframe
        plot3(or_plt, x([1:end 1], :), y([1:end 1], :), z([1:end 1], :), '-', 'Color', color_if_not_gradient, 'LineWidth', 1)
    else
        if gradient_if_not_wireframe
            axes(or_plt)
            patch(x, y, z, z, 'FaceAlpha', face_alpha)  
            colorbar(or_plt)
        else
            axes(or_plt)
            patch(x, y, z, color_if_not_gradient, 'FaceAlpha', face_alpha)
        end %if
    end %if
    plot3(xlim, [0 0], [0 0], 'k', 'LineWidth', at_origin_axis_line_widths)
    plot3([0 0], ylim, [0 0], 'k', 'LineWidth', at_origin_axis_line_widths)
    plot3([0 0], [0 0], zlim, 'k', 'LineWidth', at_origin_axis_line_widths)
end %if

if PLOT_V
    v(3*len+1:end) = [];
    cla(v_plt)
    plot(v_plt, t, v(1:3:end), 'r', t, v(2:3:end), 'g', t, v(3:3:end), 'b', 'LineWidth', LINE_WIDTH)
end %if

if PLOT_Q
    q(4*len+1:end) = [];
    cla(q_plt)
    plot(q_plt, t, q(1:4:end), 'r', t, q(2:4:end), 'g', t, q(3:4:end), 'b', t, q(4:4:end), 'k', 'LineWidth', LINE_WIDTH)
end %if

if PLOT_Q_ORB
    cla(q_orb_plt)
    plot(q_orb_plt, t, q_orb(1:4:end), 'r', t, q_orb(2:4:end), 'g', t, q_orb(3:4:end), 'b', t, q_orb(4:4:end), 'k', 'LineWidth', LINE_WIDTH)
end %if
    
if PLOT_W
    w(3*len+1:end) = [];
    cla(w_plt)
    plot(w_plt, t, w(1:3:end), 'r', t, w(2:3:end), 'g', t, w(3:3:end), 'b', 'LineWidth', LINE_WIDTH)
end %if

if PLOT_V_B
    v_b(3*len+1:end) = [];
    cla(v_b_plt)
    plot(v_b_plt, t, v_b(1:3:end), 'r', t, v_b(2:3:end), 'g', t, v_b(3:3:end), 'b', 'LineWidth', LINE_WIDTH)
end %if

if PLOT_E
    p_e(len+1:end) = [];
    cla(e_plt)
    plot(e_plt, t, p_e, 'r', 'LineWidth', LINE_WIDTH)
end %if

if PLOT_ALT
    alt(len+1:end) = [];
    cla(alt_plt)
    plot(alt_plt, t, alt, 'g', 'LineWidth', LINE_WIDTH)
end %if

if PLOT_B_STAR
    b_star(len+1:end) = [];
    cla(b_star_plt)
    plot(b_star_plt, t, b_star, 'b', 'LineWidth', LINE_WIDTH)
end %if
end %function

function fixAxes
global v_plt v_b_plt w_plt q_plt q_orb_plt e_plt alt_plt b_star_plt t PLOT_V PLOT_Q PLOT_Q_ORB PLOT_W PLOT_V_B PLOT_E PLOT_ALT PLOT_B_STAR
max_t = t(end);

if PLOT_V
    xlim(v_plt, [0 max_t])
end %if

if PLOT_Q
    xlim(q_plt, [0 max_t])
end %if

if PLOT_Q_ORB
    xlim(q_orb_plt, [0 max_t])
end %if
    
if PLOT_W
    xlim(w_plt, [0 max_t])
end %if

if PLOT_V_B
    xlim(v_b_plt, [0 max_t])
end %if

if PLOT_E
    xlim(e_plt, [0 max_t])
end %if

if PLOT_ALT
    xlim(alt_plt, [0 max_t])
end %if

if PLOT_B_STAR
    xlim(b_star_plt, [0 max_t])
end %if
end %function

function initBools
global PLOT_R PLOT_V PLOT_Q PLOT_Q_ORB PLOT_W PLOT_V_B PLOT_E PLOT_ALT PLOT_B_STAR PLOT_ORIENTATION
PLOT_R = false;
PLOT_V = false;
PLOT_Q = false;
PLOT_Q_ORB = false;
PLOT_W = false;
PLOT_V_B = false;
PLOT_E = false;
PLOT_ALT = false;
PLOT_B_STAR = false;
PLOT_ORIENTATION = false;
end %function