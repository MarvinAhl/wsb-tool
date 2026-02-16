function kep = car2keps(xx, GM)

% car2keps.m computes the Keplerian elements from the Cartesian coordinates (output 
% as a vector).
%
% PROTOTYPE:
%   kep = car2keps(xx, GM)
%
% INPUTS:
%   xx[6x1]      Position and velocity vector [km] and [km/s]
%   GM[1]        Gravitational parameter of the primary mass [km^3/s^2]
%
% OUTPUTS:
%   kep struct   Keplerian parameters in a struct "a", "e", "inc", "raan", "aop", "ta"  
%                where:
%                  - a*        Semi-major axis [km]
%                  - e         Eccentricity [-]
%                  - inc       Inclination [rad]
%                  - raan**    RAAN [rad]
%                  - aop***    Argument of pericentre [rad]
%                  - ta        True anomaly [rad]
%
%   * if e=1 (parabolic orbit) rp (pericentre radius) is calculated instead.
%
%   ** if i=0 (equatorial orbit) the node line is chosen as coincident with I.
%
%   *** if e=0 (circulat orbit) the theta=0 line is chosen as coincident with
%   the node line.
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

%% Norms of rr and vv
rr = xx(1:3);
vv = xx(4:6);

r = norm(rr);
v = norm(vv);

vr = vv'*rr/r; % radial velocity

%% Definition of h
hh = cross(rr,vv);
h = norm(hh);

%% Definition of e

ee = (1/GM)*((norm(vv)^2 - GM/r)*rr - r*vr*vv);
e = norm(ee);

%% Definition of i, OM and a

% If the orbit is parabolic a is not defined, therefore the parameter rp is
% calculated instead

if (abs(e-1) > 1e-7) %*
    E = 0.5*v^2 -GM/r;
    a = -GM/(2*E);
else
    p = h^2/GM;
    rp = p/2;
    a = rp;
end

% * A 1e-7 tollerance is used to make sure that the code doesn't break in
%case calculations introduce errors in the definition of e.

% If the orbit is equatorial (i=0) the node line is not defined. It is
% conventionally chosen as coincident with I.

i = acos(hh(3)/h);
if (abs(i) > 1e-7) %*
    N = (cross([0;0;1],hh))/norm(cross([0;0;1],hh));
    if (N(2) >= 0)
        OM = acos(N(1)); % OM is the angle between N and I
    else
        OM = 2*pi - acos(N(1));
    end
else
    N = [1;0;0];
    OM = 0;
end 

% * A 1e-7 tollerance is used to make sure that the code doesn't break in
%case calculations introduce errors in the definition of i.

%% Definition of om and theta

% If the orbit is circular (e=0) any point can be chosen as pericentre, the
% theta = 0 line is conventionally chosen as coincident with the node line.

if (abs(e) > 1e-7) %* elliptical orbit
    if (abs(i)>1e-7)
        if (ee(3) >= 0)
            om = acos(N'*ee/(e));
        else
            om = 2*pi - acos(N'*ee/(e));
        end
    else 
        % if the orbit is elliptical but equatorial ee(3)=0 so ee(2)
        % GMst be used instead.
        if (ee(2) >= 0)
            om = acos(N'*ee/(e));
        else
            om = 2*pi - acos(N'*ee/(e));
        end
    end

    edotr = clip(ee'*rr/(e*r), -1, 1);  % Sometimes larger than 1 because of numerical errors
    if (vr >= 0)
        theta = acos(edotr);
    else
        theta = 2*pi - acos(edotr);
    end
else % circular orbit
    om = 0;
    % If the orbit is circular the radial velocity is always = 0.
    % To check the angle an auxillary vector (perpendicular to N) GMst be
    % used.
    aux = cross(hh,N)/(norm(cross(hh,N)));

    if (rr'*aux >= 0)
        theta = acos(N'*rr/(r));
    else
        theta = 2*pi - acos(N'*rr/(r));
    end
end

% * A 1e-7 tollerance is used to make sure that the code doesn't break in
%case calculations introduce errors in the definition of e.

% To put angle in [0,2pi]

theta = wrapTo2Pi(theta);
om = wrapTo2Pi(om);
OM = wrapTo2Pi(OM);
i = wrapTo2Pi(i);

%% Create the struct

kep.a = a;
kep.e = e;
kep.inc = i;
kep.raan = OM;
kep.aop = om;
kep.ta = theta;

end