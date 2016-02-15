%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   KPS_Kepler2State.m
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
%   Converts Keplerian orbital elements to state vector.
%         Inputs:
%         a      -  Semi-major Axis [m]
%         e      -  Eccentricity
%         i      -  Inclination [degrees]
%         Omega  -  Longitude of Ascending Node [degrees]
%         w      -  Argument of Periapsis [degrees]
%         M      -  Mean Anomaly [degrees]
%

function [r, v] = KPS_Kepler2State(a, e, i, Omega, w, M)
if ~(isscalar(a) && isscalar(e) && isscalar(w) && isscalar(Omega) ...
        && isscalar(i) && isscalar(M))
    error('All inputs should be scalars, with angular values in degrees and a in meters.')
end %if

% Earth standard gravitational parameter
GM = 398600441800000.0;

deg_to_rad = pi/180.0;
rad_to_deg = 180.0/pi;

w = deg_to_rad * w;
Omega = deg_to_rad * Omega;
i = deg_to_rad * i;
M = deg_to_rad * M;

% numerical solve of Kepler's equation for eccentric anomaly
E_error = @(E) E - e*sin(E) - M;
E = fzero(E_error, M);

% true anomaly
nu = 2*atan2(sqrt(1+e)*sin(0.5*E), sqrt(1-e)*cos(0.5*E));

% distance to Earth center
r_c = a*(1-e*cos(E));

% orbital frame
r_o = [r_c*cos(nu), r_c*sin(nu), 0.0];
v_o = sqrt(GM*a)/r_c*[-sin(E), sqrt(1-e*e)*cos(E), 0.0];

% to ECI frame
r = [r_o(1)*(cos(w)*cos(Omega)-sin(w)*cos(i)*sin(Omega)) - r_o(2)*(sin(w)*cos(Omega)+cos(w)*cos(i)*sin(Omega))
    r_o(1)*(cos(w)*sin(Omega)+sin(w)*cos(i)*cos(Omega)) + r_o(2)*(cos(w)*cos(i)*cos(Omega)-sin(w)*sin(Omega))
    r_o(1)*sin(w)*sin(i) + r_o(2)*cos(w)*sin(i)];

v = [v_o(1)*(cos(w)*cos(Omega)-sin(w)*cos(i)*sin(Omega)) - v_o(2)*(sin(w)*cos(Omega)+cos(w)*cos(i)*sin(Omega))
    v_o(1)*(cos(w)*sin(Omega)+sin(w)*cos(i)*cos(Omega)) + v_o(2)*(cos(w)*cos(i)*cos(Omega)-sin(w)*sin(Omega))
    v_o(1)*sin(w)*sin(i) + v_o(2)*cos(w)*sin(i)];

w = w * rad_to_deg;
Omega = Omega * rad_to_deg;
i = i * rad_to_deg;
nu = nu * rad_to_deg;
M = M * rad_to_deg;
E = E * rad_to_deg;

fprintf('\nSemi-major Axis: %.14g m\n', a);
fprintf('Eccentricity: %.14g\n', e);
fprintf('Argument of Periapsis: %.14g degrees\n', w);
fprintf('Longitude of Ascending Node: %.14g degrees\n', Omega);
fprintf('Inclination: %.14g degrees\n', i);

fprintf('True Anomaly: %.14g degrees\n', nu);
fprintf('Mean Anomaly: %.14g degrees\n', M);
fprintf('Eccentric Anomaly: %.14g degrees\n', E);

fprintf('Position Vector: %.14g, %.14g, %.14g m\n', r);
fprintf('Velocity Vector: %.14g, %.14g, %.14g m/s\n\n', v);