%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   KPS_PlotPoly.m
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
%   Plots polygons of generated polygon file in 3-D to
%   examine and verify correctness.
%

function KPS_PlotPoly(polygon_file)
font_size = 24;
outer_axis_line_widths = 2;
at_origin_axis_line_widths = 1;
wireframe = false;
gradient_if_not_wireframe = true;
color_if_not_gradient = 'blue';
num_vtx = 4;
maximize_plot = true;

% for MATLAB only, won't work in Octave
face_alpha = 0.75;

if exist('OCTAVE_VERSION', 'builtin')
  face_alpha = 1.0;
  graphics_toolkit('fltk')
end %if

fid = fopen(polygon_file, 'r');
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

num_poly = size(P, 2) / num_vtx;
x = zeros(num_vtx, num_poly);
y = zeros(num_vtx, num_poly);
z = zeros(num_vtx, num_poly);
for i = 1:num_poly
    x(:, i) = P(1, num_vtx*(i-1)+1:num_vtx*i);
    y(:, i) = P(2, num_vtx*(i-1)+1:num_vtx*i);
    z(:, i) = P(3, num_vtx*(i-1)+1:num_vtx*i);
end %if
clear fid P

short_name = polygon_file;
for i = numel(polygon_file):-1:1
    if polygon_file(i) == '/' || polygon_file(i) == '\'
        short_name = polygon_file(i+1:end);
        break
    end %if
end %for

figure('Name', ['KPS - ' short_name], 'NumberTitle', 'off')
clf, hold on, grid on, axis equal
set(gca, 'FontSize', font_size, 'LineWidth', outer_axis_line_widths)
xlabel('x_b_o_d_y [m]')
ylabel('y_b_o_d_y [m]')
zlabel('z_b_o_d_y [m]')

if wireframe
    plot3(x([1:end 1], :), y([1:end 1], :), z([1:end 1], :), '-', 'Color', color_if_not_gradient, 'LineWidth', 1)
else
    if gradient_if_not_wireframe
        patch(x, y, z, z, 'FaceAlpha', face_alpha)  
        colorbar
    else
        patch(x, y, z, color_if_not_gradient, 'FaceAlpha', face_alpha)
    end %if
end %if
plot3(xlim, [0 0], [0 0], 'k', 'LineWidth', at_origin_axis_line_widths)
plot3([0 0], ylim, [0 0], 'k', 'LineWidth', at_origin_axis_line_widths)
plot3([0 0], [0 0], zlim, 'k', 'LineWidth', at_origin_axis_line_widths)
if maximize_plot
    set(gcf, 'Units', 'normalized', 'Position', [.5,0,.5,1]);
end %if