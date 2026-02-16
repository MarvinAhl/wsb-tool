function xx = kep2cars(kep, GM)

% kep2cars.m computes the Cartesian coordinates (state vector) from the Keplerian elements.
%
% PROTOTYPE:
%   xx = kep2cars(kep, GM)
%
% INPUTS:
%   kep struct with fields
%       a*[1]        Semi-major axis [km]
%       e[1]         Eccentricity [-]
%       inc[1]       Inclination [rad]
%       raan**[1]    RAAN [rad]
%       aop***[1]    Argument of pericentre [rad]
%       ta[1]        True anomaly [rad]
%   GM[1]        Gravitational parameter of the primary mass [km^3/s^2]
% 
%   * if e=1 (parabolic orbit) give rp (pericentre radius) instead.
%
% OUTPUTS:
%   xx[6x1]      Position and velocity vector [km] and [km/s]
%
% CALLED FUNCTIONS:
%   (none)
%
% CONTRIBUTORS:
%   Giacomo Burlando
%   Marta Brenna
%   Luca Matteotti
%   Marvin Ahlborn
%
% VERSIONS:
%   2025-01-07: Last version
%
% -------------------------------------------------------------------------

%% Convert variables
a = kep.a;
e = kep.e;
i = kep.inc;
OM = kep.raan;
om = kep.aop;
theta = kep.ta;

%% Creates rotation matrixes
ROM = [cos(OM), sin(OM), 0; -sin(OM), cos(OM), 0; 0, 0, 1];
Ri = [1, 0, 0; 0, cos(i), sin(i); 0, -sin(i), cos(i)];
Rom = [cos(om), sin(om), 0; -sin(om), cos(om), 0; 0, 0, 1];

T = ROM'*Ri'*Rom'; % final rotational matrix (transported)

%% Calculates p and h
if (abs(e-1) > 1e-7)
    p = a*(1-e^2);
else
    % if the orbit is parabolic, p is calculated with rp.
    rp = a;
    p = rp*2;
end
h = sqrt(p*GM);

%% Creates pericentric position and velocity vectors
rrp = p/(1+e*cos(theta))*[cos(theta); sin(theta); 0];
vvp = (GM/h)*[-sin(theta); (e+cos(theta)); 0];

%% Rotates pericentric vectors to ECI coordinates
rr = T*rrp;
vv = T*vvp;
xx = [rr; vv];

end