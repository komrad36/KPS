%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   KPS_Vis.m
%   KPS
%	
%	Author: Kareem Omar
%	kareem.omar@uah.edu
%	https://github.com/komrad36
%
%	Last updated Mar 20, 2016
%   This application is entirely my own work.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Realtime visualization of any combination of the 18 parameters
%   selectable below. Just comment out the ones you don't want.
%
%   Visualization begins in realtime mode, constantly checking
%   the outfiles for updates. When no changes occur in the files
%   for some time, KPS_Vis switches to FINAL mode, adjusting the axes
%   and allowing user interaction.
%

function KPS_Vis
% do not modify this line; configure plotting options below
[t, R, V, Q, Q_ORB, ALT, B_STAR, W, V_B, E, ORIENTATION, SEMI_MAJOR, ECC, INC, RAAN, PERIAPSIS, MEAN_ANOM, TRUE_ANOM, ECC_ANOM] = deal(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19);
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

%%% simply comment
%%% out the ones you don't want

plots = [
% R
ORIENTATION
% V
% Q
% Q_ORB
% ALT
% B_STAR
W
% V_B
% E
% SEMI_MAJOR
% ECC
% INC
% RAAN
% PERIAPSIS
% MEAN_ANOM
% TRUE_ANOM
% ECC_ANOM
];


%% Execution

if exist('OCTAVE_VERSION', 'builtin')
  face_alpha = 1.0;
  LEGEND_FONT_SIZE = 13;
  graphics_toolkit('fltk')
end %if

% number of elements for each input
elems = [1,3,3,4,4,1,1,3,3,1,0,1,1,1,1,1,1,1,1];

% names of data files
names = {'t', 'r', 'v', 'q', 'q_orb', 'alt', 'b_star', 'w', 'v_body', 'p_e'};
titles = {'', 'Position', 'Velocity', 'Quaternion to ECI Frame', 'Quaternion to Orbital Frame', 'Altitude', 'Starred Ballistic Coefficient', 'Angular Velocity', 'Velocity in Body Frame', 'Pointing Error', 'Live Orientation', 'Semi-Major Axis', 'Eccentricity', 'Inclination', 'Longitude of Ascending Node', 'Argument of Periapsis', 'Mean Anomaly', 'True Anomaly', 'Eccentric Anomaly'};
ylabels = {'', 'y_{ECI}', '[m/s]', 'Component', 'Component', '[m]', '[m^-^1]', '[s^-^1]', '[m/s]', '[\circ]', 'y_{orbital}', '[m]', [], '[\circ]', '[\circ]', '[\circ]', '[\circ]', '[\circ]', '[\circ]'};
xlabels = {'x_{ECI}', 'x_{orbital}'};
zlabels = {'z_{ECI}', 'z_{orbital}'};
legend_labels = {{}, {}, {'v_x', 'v_y', 'v_z'}, {'q_0', 'q_1', 'q_2', 'q_3'}, {'q_0', 'q_1', 'q_2', 'q_3'}, {'alt'}, {'B*'}, {'\omega_x', '\omega_y', '\omega_z'}, {'v_x', 'v_y', 'v_z'}, {'\epsilon_p'}, {}, {'a'}, {'e'}, {'i'}, {'\Omega'}, {'\omega'}, {'M'}, {'\nu'}, {'E'}};
num_plots = numel(plots);
compute_keplerian = false;

% data files required by the user's choice of plots
reqd = t;
for j = 1:num_plots
    i = plots(j);
    if i < ORIENTATION && ~any(reqd==i)
        reqd = [reqd i];
    elseif i == ORIENTATION && ~any(reqd==Q_ORB)
        reqd = [reqd Q_ORB];
        fid = fopen(poly_file, 'r');
        if fid == -1
           error('Unable to open polygon file. Aborting.')  
        end %if

        % load polygons
        % this could be textscan() in MATLAB but Octave
        % doesn't yet support '%[]' (skip characters) arguments
        % so it's handrolled for now
        P = zeros(3, 0);
        k = 1;
        while true
          ln = fgetl(fid);
          if ln == -1, break, end

          if numel(ln) == 0 || ln(1) == '#', continue, end
          tokens = str2double(strsplit(ln, ','));
          if numel(tokens) ~= 3
            error('Error parsing file!')
          end
          P(:, k) = tokens;
          k = k + 1;
        end %while

        fclose(fid);
        clear fid
    else
        if ~any(reqd==R), reqd = [reqd R]; end
        if ~any(reqd==V), reqd = [reqd V]; end
        compute_keplerian = true;
    end %if
end %for

num_data = numel(reqd);
h_f = zeros(1, 19);
data = cell(1, 19);
for j = 1:num_data
    i = reqd(j);
    h_f(i) = fopen(strcat(names{i}, '.bin'), 'r');
    if h_f(i) == -1, error('Failed to open file.'), end
end %for

h_fig = figure('Name', 'KPS - REALTIME', 'NumberTitle', 'off'); clf
if MAXIMIZE_PLOT, set(h_fig, 'units','normalized','outerposition',[0 0 1 1]), end

f_blank = {0
        @() plot3([0 0; 0 0], [0 0; 0 0], [0 0; 0 0], 'g', 'LineWidth', 2)
        @() plot(0, 0, 'r', 0, 0, 'g', 0, 0, 'b', 'LineWidth', LINE_WIDTH)
        @() plot(0, 0, 'r', 0, 0, 'g', 0, 0, 'b', 0, 0, 'k', 'LineWidth', LINE_WIDTH)
        @() plot(0, 0, 'r', 0, 0, 'g', 0, 0, 'b', 0, 0, 'k', 'LineWidth', LINE_WIDTH)
        @() plot(0, 0, 'g', 'LineWidth', LINE_WIDTH)
        @() plot(0, 0, 'b', 'LineWidth', LINE_WIDTH)
        @() plot(0, 0, 'r', 0, 0, 'g', 0, 0, 'b', 'LineWidth', LINE_WIDTH)
        @() plot(0, 0, 'r', 0, 0, 'g', 0, 0, 'b', 'LineWidth', LINE_WIDTH)
        @() plot(0, 0, 'r', 'LineWidth', LINE_WIDTH)
        []
        @() plot(0, 0, 'g', 'LineWidth', LINE_WIDTH)
        @() plot(0, 0, 'b', 'LineWidth', LINE_WIDTH)
        @() plot(0, 0, 'r', 'LineWidth', LINE_WIDTH)
        @() plot(0, 0, 'g', 'LineWidth', LINE_WIDTH)
        @() plot(0, 0, 'b', 'LineWidth', LINE_WIDTH)
        @() plot(0, 0, 'r', 'LineWidth', LINE_WIDTH)
        @() plot(0, 0, 'g', 'LineWidth', LINE_WIDTH)
        @() plot(0, 0, 'b', 'LineWidth', LINE_WIDTH)
};

% init plots
h_plots = zeros(1, 11);
h_lines = cell(1, 11);
for j = 1:num_plots
    i = plots(j);
    h_plots(i) = subplot(num_plots, 1, j);
    h_lines{i} = f_blank{i}();
    if numel(legend_labels{i})
        h_legend = legend(legend_labels{i}, 'location', 'eastoutside');
        % font size done separately due to strange Octave bug when combining
        % cell arrays and multiple parameters
        set(h_legend, 'FontSize', LEGEND_FONT_SIZE)
    end %if
    if i == R
        box off
        grid on
        R_str = num2str(R);
        set(h_lines{i}(1), 'xdatasource', ['data{' R_str '}(1:3:end)'], 'ydatasource', ['data{' R_str '}(2:3:end)'], 'zdatasource', ['data{' R_str '}(3:3:end)'])
        set(h_lines{i}(2), 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'MarkerSize', 50);
    end %if
    if i ~= R && i ~= ORIENTATION
        for k = 1:numel(h_lines{i})
          set(h_lines{i}(k), 'xdatasource', 'data{1}', 'ydatasource', ['data{' num2str(i) '}(' num2str(k) ':' num2str(elems(i)) ':end)'])
        end %for
    end %if
    title(titles{i}, 'FontSize', FONT_SIZE)
    ylabel(ylabels{i}, 'FontSize', FONT_SIZE)
    set(gca, 'FontSize', FONT_SIZE, 'LineWidth', AXIS_LINE_WIDTH);
    if i == E, ylim([0 180]), end
    if i == INC, ylim([-90 90]), end
    if i >= RAAN, ylim([0 360]), end
    if i == R || i == ORIENTATION
        axis equal
        xlabel(xlabels{1+(i~=R)})
        zlabel(zlabels{1+(i~=R)})
        if i == ORIENTATION
            hold on
            colorbar
        end %if
    end %if
end %for

% print xlabel if there are plots other than the 3-D ones, orientation and position 
if ~((num_plots == 1 && (plots(1) == ORIENTATION || plots(1) == R)) || (num_plots == 2 && any(plots == ORIENTATION) && any(plots == R)))
    xlabel('Time [s]', 'FontSize', FONT_SIZE)
end %if
h_text = text(mean(xlim), mean(ylim), 'Loading...', 'FontSize', 50, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');

ticks_without_change = 0;
rad_to_deg = 180.0/pi;
while true
    if ~ishghandle(h_fig), break, end
    drawnow
    
    old_t_size = numel(data{1});
    % update data
    for j = 1:num_data
        i = reqd(j);
        fseek(h_f(i), 8*numel(data{i}), 'bof');
        data{i} = [data{i}; fread(h_f(i), Inf, 'double')];
    end %for
    
    % shorten all vectors to length of shortest one
    % minus 1, since with realtime data the very last
    % entry might still have been incomplete and thus
    % corrupt
    min_len = floor(min(bsxfun(@rdivide, cellfun(@numel, data(reqd)), elems(reqd))))-1;
    for j = 1:num_data
        i = reqd(j);
        data{i}(elems(i)*min_len+1:end) = [];
    end %for
       
    if min_len == old_t_size
        ticks_without_change = ticks_without_change + 1;
        if ticks_without_change > 20
            % no change for a while
            % switch to Final mode
            if ~ishghandle(h_fig), break, end
            for j = 1:num_plots
                i = plots(j);
                if i ~= R && i ~= ORIENTATION, xlim(h_plots(i), [0 data{1}(end)]), end
            end %for
            set(h_fig, 'Name', 'KPS - Final')
            break
        end %if
    else
        ticks_without_change = 0;
        
        % do Keplerian if necessary
        if compute_keplerian
            rng_a = numel(data{SEMI_MAJOR})+1;
            r = reshape(data{R}(3*rng_a-2:end), 3, []);
            v = reshape(data{V}(3*rng_a-2:end), 3, []);
            r_mag = sqrt(sum(r.*r));
            v_mag = sqrt(sum(v.*v));

            % orbital momentum
            h = cross(r, v);

            % eccentricity vector
            % the magic num is 1/mu, Earth's standard gravitational parameter,
            % precomputed for speed
            e_vec = 2.508777951886354381858e-15*cross(v, h)-bsxfun(@rdivide, r, r_mag);
            % scalar eccentricity
            e = sqrt(sum(e_vec.*e_vec));

            % ascending node vector
            n_1 = -h(2, :);
            n_2 = h(1, :);
            % n_3 is 0
            n_mag = hypot(n_1, n_2);

            % true anomaly
            nu = real(acos(dot(e_vec, r)./(e.*r_mag)));
            idx = dot(r, v) < 0;
            nu(idx) = 2*pi - nu(idx);

            % inclination
            inc = real(acos(h(3, :)./sqrt(sum(h.*h))));

            % eccentric anomaly
            E = real(2*atan2(tan(0.5*nu), sqrt((1+e)./(1-e))));
            idx = E < 0;
            E(idx) = E(idx) + 2*pi;

            % RAAN
            Omega = real(acos(n_1./n_mag));
            idx = n_2 < 0;
            Omega(idx) = 2*pi - Omega(idx);

            % argument of periapsis
            w = real(acos((n_1.*e_vec(1,:) + n_2.*e_vec(2,:))./(n_mag.*e)));
            idx = e_vec(3, :) < 0;
            w(idx) = 2*pi - w(idx);

            % mean anomaly
            M = E - e.*sin(E);
            idx = M < 0;
            M(idx) = M(idx) + 2*pi;

            rng_b = rng_a + size(r, 2) - 1;
            % the magic num is 1/mu, Earth's standard gravitational parameter,
            % precomputed for speed
            data{SEMI_MAJOR}(rng_a:rng_b) = 1./(2./r_mag-v_mag.*2.508777951886354381858e-15.*v_mag);
            data{TRUE_ANOM}(rng_a:rng_b) = nu*rad_to_deg;
            data{INC}(rng_a:rng_b) = inc*rad_to_deg;
            data{ECC}(rng_a:rng_b) = e;
            data{RAAN}(rng_a:rng_b) = Omega*rad_to_deg;
            data{PERIAPSIS}(rng_a:rng_b) = w*rad_to_deg;
            data{MEAN_ANOM}(rng_a:rng_b) = M*rad_to_deg;
            data{ECC_ANOM}(rng_a:rng_b) = E*rad_to_deg;
        end %if
        
        if ishandle(h_text), delete(h_text), end
    end %if
   
    if ~ishghandle(h_fig), break, end
    
    % plot updated data
    for j = 1:num_plots
        i = plots(j);
        if i == ORIENTATION
            P_rot = zeros(size(P));
            for k = 1:size(P, 2)
                P_rot(:, k) = P(:, k) + cross(2*data{Q_ORB}(end-2:end), cross(data{Q_ORB}(end-2:end), P(:, k))+data{Q_ORB}(end-3)*P(:, k));
            end %for
            num_poly = size(P, 2) / num_vtx;
            x = zeros(num_vtx, num_poly);
            y = zeros(num_vtx, num_poly);
            z = zeros(num_vtx, num_poly);
            for k = 1:num_poly
                x(:, k) = P_rot(1, num_vtx*(k-1)+1:num_vtx*k);
                y(:, k) = P_rot(2, num_vtx*(k-1)+1:num_vtx*k);
                z(:, k) = P_rot(3, num_vtx*(k-1)+1:num_vtx*k);
            end %if
            cla(h_plots(i))
            if wireframe
                plot3(h_plots(i), x([1:end 1], :), y([1:end 1], :), z([1:end 1], :), '-', 'Color', color_if_not_gradient, 'LineWidth', 1)
            else
                if gradient_if_not_wireframe
                    patch(x, y, z, x, 'FaceAlpha', face_alpha, 'Parent', h_plots(i))  
                else
                    patch(x, y, z, color_if_not_gradient, 'FaceAlpha', face_alpha, 'Parent', h_plots(i))
                end %if
            end %if
            plot3(h_plots(i), xlim(h_plots(i)), [0 0], [0 0], 'k', 'LineWidth', at_origin_axis_line_widths)
            plot3(h_plots(i), [0 0], ylim(h_plots(i)), [0 0], 'k', 'LineWidth', at_origin_axis_line_widths)
            plot3(h_plots(i), [0 0], [0 0], zlim(h_plots(i)), 'k', 'LineWidth', at_origin_axis_line_widths)
        end %if
        if ~exist('OCTAVE_VERSION', 'builtin')
            refreshdata(h_plots(i), 'caller')
        end %if
    end %for
    if exist('OCTAVE_VERSION', 'builtin')
        refreshdata(h_fig, 'caller')
    end %if
end %while
arrayfun(@fclose, h_f(reqd));
end %function