%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   KPS_GenPoly.m
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
%   Generates a polygon configuration file for input to the KPS
%   propagation framework. The idea is that users can modify this script,
%   or copies of it, to be arbitrarily complex, allowing easy modification
%   of a parameter such as a solar panel deployment angle without having to
%   manually recompute numerical values for all vertices.
%
%   For example, this script currently writes a 6U "dart" configuration
%   satellite, and rather than having to manually compute vertices based
%   on the panel length and angle, simply modify the 'H' and 'PA' variables
%   and run this script to generate the updated polygon file.
%   
%   Verify correct polygon generation with KPS_PlotPoly.
%

function KPS_GenPoly(outfile)

num_vtx = 4;

% folding panel length [m]
H = 0.3;

% panel angle
PA = 0.0 / 180.0 * pi;

% header string to go at the top of the output polygon file
header = ['# 6U Dart\n# Panel Length: ', num2str(H), ...
    '\n# Panel Angle: ', num2str(180.0/pi*PA)];

% The actual body frame satellite as a series of polygons.
% This is a 6U 'dart' configuration facing in the +z with all 4 panels
% bilaterally symmetrically deployed to the panel angle specified above.
poly = [
    
0.1, -0.05, -0.15
0.1, 0.05, -0.15
0.1, 0.05, 0.15
0.1, -0.05, 0.15
-0.1, -0.05, -0.15
-0.1, 0.05, -0.15
-0.1, 0.05, 0.15
-0.1, -0.05, 0.15
-0.1, 0.05, -0.15
0.1, 0.05, -0.15
0.1, 0.05, 0.15
-0.1, 0.05, 0.15
-0.1, -0.05, -0.15
0.1, -0.05, -0.15
0.1, -0.05, 0.15
-0.1, -0.05, 0.15
-0.1, -0.05, 0.15
0.1, -0.05, 0.15
0.1, 0.05, 0.15
-0.1, 0.05, 0.15
-0.1, -0.05, -0.15
0.1, -0.05, -0.15
0.1, 0.05, -0.15
-0.1, 0.05, -0.15
0.1, 0.05, -0.15
0.1, -0.05, -0.15
0.1 + H * sin(PA), -0.05, -0.15 + H * cos(PA)
0.1 + H * sin(PA), 0.05, -0.15 + H * cos(PA)
-0.1, 0.05, -0.15
-0.1, -0.05, -0.15
-0.1 - H * sin(PA), -0.05, -0.15 + H * cos(PA)
-0.1 - H * sin(PA), 0.05, -0.15 + H * cos(PA)
0.1, 0.05, -0.15
-0.1, 0.05, -0.15
-0.1, 0.05 + H * sin(PA), -0.15 + H * cos(PA)
0.1, 0.05 + H * sin(PA), -0.15 + H * cos(PA)
0.1, -0.05, -0.15
-0.1, -0.05, -0.15
-0.1, -0.05 - H * sin(PA), -0.15 + H * cos(PA)
0.1, -0.05 - H * sin(PA), -0.15 + H * cos(PA)

];

fid = fopen(outfile, 'w');
if (fid == -1)
    error('Failed to open file.')
end %if

% write polygons to file
fprintf(fid, strcat(header, '\n\n'));
for i = 1:num_vtx:size(poly, 1)
    fprintf(fid, '%.16g, %.16g, %.16g\n', poly(i:i+num_vtx-1, :)');
    fprintf(fid, '\n');
end %for
fclose(fid);

fprintf('Successfully wrote %s.\n', outfile)