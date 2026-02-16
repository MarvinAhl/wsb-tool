% phi is the angle of the earth-moon line wrt to the vernal equinox
% th is the angle of the sun wrt the negative earth-moon line (measured clockwise)

function [phi, th] = epoch_to_BCRFBP_angles(et)
    frame = 'ECLIPJ2000';
    sun_label = 'Sun';
    earth_label = 'Earth';
    moon_label = 'Moon';

    rr_sun = cspice_spkpos(sun_label, et, frame, 'NONE', earth_label);
    rr_moon = cspice_spkpos(moon_label, et, frame, 'NONE', earth_label);

    phi = atan2(rr_moon(2), rr_moon(1));

    sun_ang = atan2(rr_sun(2), rr_sun(1));
    th = wrapTo2Pi(phi - sun_ang + pi);
end