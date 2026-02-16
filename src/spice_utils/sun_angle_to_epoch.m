% Simple iterative scheme to find closest epoch that has sun angle th_target
function et = sun_angle_to_epoch(th_target, et0, mdata)
    time_per_angle = 1 / (mdata.omega_s * mdata.uangfreq_sinv);

    last_dth = inf;
    next_et = et0;

    tol = deg2rad(0.01);
    while abs(last_dth) > tol
        [~, th] = epoch_to_BCRFBP_angles(next_et);
        dth = wrapToPi(th - th_target);

        next_et = next_et - time_per_angle * dth;
        last_dth = dth;
    end

    et = next_et;
end