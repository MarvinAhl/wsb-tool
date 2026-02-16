% Transforms from inertial ECLIPJ2000 to either inertial J2000 or earth
% fixed IAU_EARTH

function traj_ece = traj_ecl_to_equ(ets, traj_eci, earth_fixed)
    to_frame = 'J2000';
    if earth_fixed
        to_frame = 'ITRF93';
    end

    N_s_els = size(traj_eci, 2);

    Ts_ECIN_ECEF = [];
    if N_s_els == 3
        Ts_ECIN_ECEF = cspice_pxform('ECLIPJ2000', to_frame, ets);
    elseif N_s_els == 6
        Ts_ECIN_ECEF = cspice_sxform('ECLIPJ2000', to_frame, ets);
    end

    traj_ece = zeros(length(ets), N_s_els);
    for ii = 1 : length(ets)
        T = Ts_ECIN_ECEF(:, :, ii);
        traj_ece(ii, :) = T * traj_eci(ii, :)';
    end
end