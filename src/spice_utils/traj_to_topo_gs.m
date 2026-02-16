function [traj_rho, traj_az, traj_el] = traj_to_topo_gs(ets, traj_ECLJ2000, gs_code)
    if size(ets, 1) > size(ets, 2)
        ets = ets';
    end

    gs_name = strcat('NDOSL_', gs_code);
    gs_frame = strcat('NDOSL_', gs_code, '_TOPO');

    % Change frame to ground station topocentric
    Ts_ECLJ2000_gs = cspice_pxform('ECLIPJ2000', gs_frame, ets);

    traj_gs = zeros(size(traj_ECLJ2000, 1), 3);
    for ii = 1 : length(ets)
        T = Ts_ECLJ2000_gs(:, :, ii);
        traj_gs(ii, :) = T * traj_ECLJ2000(ii, 1:3)';
    end

    % Change origin
    gs_pos = cspice_spkpos(gs_name, ets(1), gs_frame, 'NONE', 'Earth');  % Only ets(1) because n gs_frame, gs_name is time independent
    traj_gs = traj_gs - gs_pos';

    % Convert to latitudinal coordinates
    [traj_rho, traj_az, traj_el] = cspice_reclat(traj_gs');
end