%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   KPS_GenOrbit.m
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
%   Generates a satellite initial position and velocity
%   in the format required by KPS such that the satellite
%   is in an x by y kilometer orbit at i degrees inclination, where
%   the user specifies x, y, and i. The satellite begins at the 'x'
%   portion of the orbit.
%

function KPS_GenOrbit(x, y, inc)

% Earth standard gravitational parameter
GM = 398600441800000.0;

% Earth radius [m]
R  = 6371000.0;

M_PER_KM = 1000.0;

min_alt = 100;
max_alt = 1000;

if ~isscalar(x) || ~isscalar(y) || ~isscalar(inc)
    error('Please enter scalar values only.')
end %if

if x < min_alt || y < min_alt
    error(['Altitudes must be ' num2str(min_alt) ' km or above. Aborting.'])
end %if

if x > max_alt || y > max_alt
    fprintf(['\nWARNING: KPS is intended for altitudes below ' ...
        num2str(max_alt) ' km.\n' ...
        'Proceed at your own risk.\n'])
end %if

if inc < -90.0 || inc > 90.0
    error(['Inclination must be -90 to 90' deg_symbol '. Aborting.'])
end %if

x_r = R + M_PER_KM*x;
y_r = R + M_PER_KM*y;

init_pos = [x_r, 0.0, 0.0];

% semi-major axis
a = 0.5*(x_r + y_r);

% vis-viva
init_vel_mag = sqrt(GM*(2/x_r-1/a));
init_vel = [0.0, init_vel_mag*cosd(inc), init_vel_mag*sind(inc)];

fprintf('\nGenerating stats for a %.14gx%.14g km orbit at %.14g degrees...\n', ...
    x, y, inc)
fprintf('\nSAT_INIT_POS = %.14g, %.14g, %.14g\n', init_pos)
fprintf('SAT_INIT_V = %.14g, %.14g, %.14g\n\n', init_vel)