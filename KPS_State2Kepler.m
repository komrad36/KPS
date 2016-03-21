%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   KPS_State2Kepler.m
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
%   Converts state vectors to Keplerian orbital elements.
%
%   A state vector may be directly supplied, OR a single number
%   can be supplied - a time in seconds. KPS_State2Kepler will
%   search through the most recent KPS run and retrieve the
%   first state vector AFTER that time for conversion.
%

function [a, e, w, Omega, i, M] = KPS_State2Kepler(r, v)
% if in time retrieval mode
if isscalar(r) && nargin == 1
    % grab data from files
    t_val = r;
    f_t = fopen('t.bin', 'r');
    if (f_t == -1)
        error('Failed to open file.')
    end %if
    t = fread(f_t, Inf, 'double');
    
    % discard last entry in case it is corrupt
    % (if it's being written realtime)
    t = t(1:end-1);

    f_r = fopen('r.bin', 'r');
    if (f_r == -1)
        error('Failed to open file.')
    end %if
    r = fread(f_r, Inf, 'double');
    r = r(1:3*ceil(numel(r)/3) - 3);

    f_v = fopen('v.bin', 'r');
    if (f_v == -1)
        error('Failed to open file.')
    end %if
    v = fread(f_v, Inf, 'double');
    v = v(1:3*ceil(numel(v)/3) - 3);
    
    fclose('all');
    
    if t_val >= t(1)
        for i = 1:numel(t)
            if t(i) >= t_val
                % found the first time step after requested time
                r = r(3*i - 2:3*i);
                v = v(3*i - 2:3*i);
                fprintf('\nTime: %.14g s\n', t(i))
                clear t
                break
            end %if
        end %for
    end %if
    if exist('t', 'var')
        error('Requested time not in range of files!')
    end %if
end %if
    
if numel(r) ~= 3 || numel(v) ~= 3
    error('r and v should be Cartesian orbital state 3-vectors, units of meters')
end %if

% Earth gravitational parameter
GM = 398600441800000.0;

rad_to_deg = 180.0/pi;

v_mag = norm(v);
r_mag = norm(r);

% orbital momentum
h = cross(r, v);

% eccentricity vector
e_vec = cross(v, h)/GM-r/r_mag;

% scalar eccentricity
e = norm(e_vec);

% ascending node vector
n = [-h(2), h(1), 0.0];
n_mag = norm(n);

% true anomaly
nu = real(acos(dot(e_vec, r)/(e*r_mag)));
if dot(r, v) < 0.0
    nu = 2*pi - nu;
end %if

% inclination
i = real(acos(h(3)/norm(h)));

% eccentric anomaly
E = real(2*atan2(tan(0.5*nu), sqrt((1+e)/(1-e))));
if E < 0
    E = E + 2*pi;
end %if

% RAAN
Omega = real(acos(n(1)/n_mag));
if n(2) < 0.0
    Omega = 2*pi - Omega;
end %if

% argument of periapsis
w = real(acos(dot(n, e_vec)/(n_mag*e)));
if e_vec(3) < 0.0
    w = 2*pi - w;
end %if

% mean anomaly
M = E - e*sin(E);
if M < 0
    M = M + 2*pi;
end %if

% semi-major axis
a = 1.0/(2.0/r_mag-v_mag*v_mag/GM);

nu = nu*rad_to_deg;
i = i*rad_to_deg;
E = E*rad_to_deg;
Omega = Omega*rad_to_deg;
w = w*rad_to_deg;
M = M*rad_to_deg;

fprintf('\nPosition Vector: %.14g, %.14g, %.14g m\n', r)
fprintf('Velocity Vector: %.14g, %.14g, %.14g m/s\n', v)

fprintf('True Anomaly: %.14g degrees\n', nu)
fprintf('Mean Anomaly: %.14g degrees\n', M)
fprintf('Eccentric Anomaly: %.14g degrees\n', E)

fprintf('Semi-major Axis: %.14g m\n', a)
fprintf('Eccentricity: %.14g\n', e)
fprintf('Argument of Periapsis: %.14g degrees\n', w)
fprintf('Longitude of Ascending Node: %.14g degrees\n', Omega)
fprintf('Inclination: %.14g degrees\n\n', i)
